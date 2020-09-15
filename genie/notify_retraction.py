import logging

import pandas as pd
import synapseclient
from synapseclient import Synapse

from . import dashboard_table_updater

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_latest_public_release(syn: Synapse, release: str) -> str:
    """Finds latest public release

    Args:
        syn: Synapse connection
        release: GENIE release version

    Returns:
        Synapse folder id

    """
    major_release = release.split(".")[0]
    public_major_release = int(major_release) - 1
    public_rel = syn.tableQuery(
        "select * from syn22233011 where "
        f"name like 'Release {public_major_release}.%-public'"
    )
    public_releasedf = public_rel.asDataFrame()
    if public_releasedf.empty:
        raise ValueError("Could not find public release "
                         f"for: {public_major_release}")
    return public_releasedf['id'][0]


def find_release(syn: Synapse, release: str) -> str:
    """Finds the Synapse id of a private consortium release folder

    Args:
        syn: Synapse connection
        release: GENIE release version

    Returns:
        Synapse folder id

    """
    release_synid = syn.tableQuery(
        "select distinct(parentId) from syn16804261 where "
        f"release = '{release}'"
    )
    releasedf = release_synid.asDataFrame()
    if releasedf.empty:
        raise ValueError("Please specify correct release value")
    return releasedf.iloc[0,0]


def notify(syn: Synapse, possible_retracted: pd.DataFrame):
    """Notify centers of samples to be retracted

    Args:
        syn: Synapse connection
        possible_retracted: Possibly retracted samples dataframe

    """
    possible_retracted['CENTER'] = [
        assay.split("-")[0] for assay in possible_retracted['SEQ_ASSAY_ID']
    ]
    logger.info(possible_retracted['CENTER'].value_counts())

    center_retracted = possible_retracted.groupby("CENTER")

    center_participants = syn.tableQuery(
        "select username, center from syn18344792"
    )
    center_participantsdf = center_participants.asDataFrame()

    for center, df in center_retracted:
        implicit_retracted = df['SAMPLE_ID'][df['RETRACTION'] == "implicit"]
        implicit_sample_str = ",".join(implicit_retracted)

        explicit_retracted = df['SAMPLE_ID'][df['RETRACTION'] == "explicit"]
        explicit_sample_str = ",".join(explicit_retracted)

        retraction_email = (
            f"Dear {center},\n\n"
            "This is to notify you that your samples will be retracted "
            "from the latest public release! "
            "Please respond to this email to confirm if the below samples "
            "should be retracted!\n\n"
            f"{len(implicit_retracted)} implicitly retracted samples:\n"
            f"{implicit_sample_str}\n\n"
            f"{len(explicit_retracted)} explicitly retracted samples:\n"
            f"{explicit_sample_str}\n\n"
            "Best,\nGENIE administrator"
        )
        center_users_idx = center_participantsdf['center'] == center
        center_users = center_participantsdf['username'][center_users_idx]
        
        emaillist = [3324230, 1968150]
        emaillist.extend(center_users.tolist())
        # Remove this when finalized
        emaillist = [3324230]
        syn.sendMessage(userIds=emaillist,
                        messageSubject=f"GENIE retraction policy: {center}",
                        messageBody=retraction_email)


def annotate_with_retraction_type(
        syn: Synapse,
        possible_retracted: pd.DataFrame
    ) -> pd.DataFrame:
    """Determine if sample is explcitly or implicitly retracted

    Args:
        syn: Synapse connection
        possible_retracted: Possibly retracted samples dataframe

    Returns:
        Same dataframe with retraction type column

    """
    # assume most samples are implicitly retracted
    possible_retracted['RETRACTION'] = "implicit"
    # Get explicitly retracted
    retracted_patient = syn.tableQuery(
        "select geniePatientId from syn11564409"
    )
    retracted_patientdf = retracted_patient.asDataFrame()
    retracted_sample = syn.tableQuery("select genieSampleId from syn8534758")
    retracted_sampledf = retracted_sample.asDataFrame()

    retracted_patients_idx = possible_retracted['PATIENT_ID'].isin(
        retracted_patientdf['geniePatientId']
    )
    retracted_samples_idx = possible_retracted['SAMPLE_ID'].isin(
        retracted_sampledf['genieSampleId']
    )

    possible_retracted['RETRACTION'][
        retracted_patients_idx & retracted_samples_idx
    ] = 'explicit'
    return possible_retracted


def get_possible_retracted(syn: Synapse,
                           public_clindf: pd.DataFrame,
                           release_clindf: pd.DataFrame) -> pd.DataFrame:
    """Get samples that are possibly retracted
    The reason 'possibly' is because centers can rename their samples

    Args:
        syn: Synapse connection
        public_clindf: Public release clinical dataframe
        release_clindf: Specified release clinical dataframe

    Returns:
        Possibly retracted samples dataframe

    """
    # Only samples that are not in the database any longer are the ones that
    # are excluded
    sample_table = syn.tableQuery("SELECT SAMPLE_ID FROM syn7517674")
    db_samplesdf = sample_table.asDataFrame()

    # check for samples in public release that are no longer in
    # consortium release
    exist_idx = public_clindf['SAMPLE_ID'].isin(release_clindf['SAMPLE_ID'])
    in_db_samples_idx = public_clindf['SAMPLE_ID'].isin(
        db_samplesdf['SAMPLE_ID']
    )
    # Get all possibly retracted
    # The reason possibly because samples could have been renamed
    possible_retracted = public_clindf[~exist_idx & ~in_db_samples_idx]

    return possible_retracted


def main(syn: Synapse, release: str):
    """Notify centers of samples to be retracted

    Args:
        syn: Synapse connect
        release: GENIE release version

    """
    # Get latest public release folder synapse id and clinical file
    pub_release_synid = get_latest_public_release(syn, release)
    # Get latest release synapse id
    release_synid = find_release(syn, release)

    public_file_mapping = dashboard_table_updater.get_file_mapping(
        syn, pub_release_synid
    )
    release_file_mapping = dashboard_table_updater.get_file_mapping(
        syn, release_synid
    )
    # Read in clinical files
    public_clin_ent = syn.get(public_file_mapping['clinical'])
    public_clindf = pd.read_csv(public_clin_ent.path, sep="\t", comment="#")
    release_clin_ent = syn.get(release_file_mapping['clinical'],
                               followLink=True)
    release_clindf = pd.read_csv(release_clin_ent.path, sep="\t", comment="#")
    # Get retracted samples
    possible_retracted = get_possible_retracted(syn, public_clindf,
                                                release_clindf)
    # Get retraction type of samples
    possible_retracted = annotate_with_retraction_type(syn, 
                                                       possible_retracted)
    # Notify sites to confirm retraction
    notify(syn, possible_retracted)
