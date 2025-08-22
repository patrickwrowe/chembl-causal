import pandas as pd
import numpy as np
import sqlite3

CHEMBLDB_PATH = "/Users/patrickrowe/Documents/Code/data/chembl/chembl_35/chembl_35_sqlite/chembl_35.db"

# Incomplete remapping of diverse chembl assay descriptions to simplified versions.
# Many to one remapping is used to combine (probably) the same assay values into single cols
# DON'T use for production without checking, none of this has been carefully verified.
CHEMBL_ASSAY_DESCRIPTION_REMAPPING = {
    # pKa assays
    'value_ASTRAZENECA: Most basic pKa value (pKa B1) determined by absorption and potentiometric titration using standard methodology from Sirius Analytical. Experimental range Bases: >= 2.': 'pKa',
    'value_Dissociation constant (pKa)': 'pKa',
    'value_Dissociation constant, pKa of the compound': 'pKa',
    'value_Ionization constant (pKa)': 'pKa',
    
    # LogD (distribution coefficient at pH 7.4)
    'value_ASTRAZENECA: Octan-1-ol/water (pH7.4) distribution coefficent measured  by a shake flask method described in J. Biomol. Screen. 2011, 16, 348-355. Experimental range -1.5 to 4.5': 'LogD',
    'value_Distribution coefficient, log D of the compound': 'LogD',
    'value_Distribution coefficient, log D of the compound at pH 7.4': 'LogD',
    'value_Lipophilicity, log D at pH 7.4': 'LogD',
    'value_Lipophilicity, log D of compound at pH 7.4': 'LogD',
    'value_Lipophilicity, log D of the compound': 'LogD',
    'value_Lipophilicity, log D of the compound at pH 7.4': 'LogD',
    'value_Lipophilicity, log D of the compound at pH 7.4 by HPLC analysis': 'LogD',
    'value_Lipophilicity, log D of the compound at pH 7.4 by chromatographic method': 'LogD',
    'value_Lipophilicity, log D of the compound at pH 7.4 by shake flask method': 'LogD',
    'value_Lipophilicity, logD of compound at pH 7.4': 'LogD',
    'value_Lipophilicity, logD of the compound at pH 7.4': 'LogD',
    'value_Partition coefficient (logD7.4)': 'LogD',
    'value_SUPPLEMENTARY: Lipophilicity (log D) Determination by high-throughput shake-flask': 'LogD',
    
    # LogP (partition coefficient)
    'value_Calculated partition coefficient (clogP)': 'cLogP',
    'value_Lipophilicity, log P of the compound': 'LogP',
    'value_Lipophilicity, logP of the compound': 'LogP',
    'value_Octanol-buffer partition coefficient, log P of the compound': 'LogP',
    'value_Octanol-water partition coefficient, log P of the compound': 'LogP',
    'value_Partition coefficient (logP)': 'LogP',
    'value_Partition coefficient (logP) (HPLC)': 'LogP',
    'value_Partition coefficient, log P of the compound': 'LogP',
    
    # Aqueous solubility (pH 7.4)
    'value_ASTRAZENECA: Solubility in pH7.4 buffer using solid starting material using the method described in J. Assoc. Lab. Autom. 2011, 16, 276-284. Experimental range 0.10 to 1500 uM': 'Aqueous Solubility',
    'value_Aqueous solubility': 'Aqueous Solubility',
    'value_Aqueous solubility at pH 7.4': 'Aqueous Solubility',
    'value_Aqueous solubility of compound': 'Aqueous Solubility', # Assuming this is in water....,
    'value_Aqueous solubility of the compound': 'Aqueous Solubility', # Water?
    'value_Aqueous solubility of the compound at pH 7.4': 'Aqueous Solubility',
    'value_Solubility at pH 7.4': 'Aqueous Solubility',
    'value_Solubility in water': 'Aqueous Solubility',
    'value_Solubility of compound in water': 'Aqueous Solubility',
    'value_Solubility of the compound': 'Aqueous Solubility',
    'value_Solubility of the compound in water': 'Aqueous Solubility',
    'value_Solubility of the compound at pH 7.4': 'Aqueous Solubility',
    
    # Kinetic solubility
    'value_Kinetic aqueous solubility of the compound': 'Kinetic Solubility',
    'value_Kinetic solubility of compound': 'Kinetic Solubility',
    'value_Kinetic solubility of the compound': 'Kinetic Solubility',
    'value_Kinetic solubility of the compound at pH 7.4': 'Kinetic Solubility',
    
    # HPLC retention time
    'value_Retention time of the compound': 'Retention Time',
    'value_Retention time of the compound by HPLC analysis': 'Retention Time',
    
    # Special pH conditions
    'value_Partition coefficient (logD6.5)': 'LogD_pH6.5',
    'value_Solubility of the compound at pH 6.8': 'Solubility_pH6.8',
    'value_Solubility of the compound at pH 7': 'Solubility_pH7',
    
    # DMSO solubility
    'value_Solubility of the compound in DMSO': 'DMSO Solubility',
    
    # Specialized assays
    'value_SUPPLEMENTARY: Lyophilisation Solubility Assay (LYSA)': 'LYSA',
    'value_SUPPLEMENTARY: PAMPA permeability assay': 'PAMPA Permeability'
}

def chembl_physicochemical_assays(min_raw_assay_measurements: int = 500, min_measurements_per_mol: int = 2, atol: float = 1e-9, verbose: bool = False):

    SQL_QUERY_WITH_IDS = """SELECT 
                        molregno, assays.assay_id, description, confidence_score, chembl_id, data_validity_comment, value, units, type
                        FROM assays
                        INNER JOIN activities AS act
                        ON assays.assay_id = act.assay_id 
                        WHERE ASSAY_TYPE is 'P'
                    """
    
    # do dat query
    con = sqlite3.connect(CHEMBLDB_PATH)
    chembl_a_p_raw = pd.read_sql_query(SQL_QUERY_WITH_IDS, con=con)
    con.close()

    # Only select assays with 
    chembl_a_p = chembl_a_p_raw[chembl_a_p_raw.groupby("description").transform("size") > min_raw_assay_measurements]

    pivoted = pd.pivot_table(chembl_a_p, values=['value'], index='molregno', columns=['description']).sort_index().pipe(
        lambda s: s.set_axis(s.columns.map("_".join), axis=1)
    )

    # Remap columns which _probably_ refer to identical assays, and take the mean where multiple values are found
    pivoted = pivoted.rename(columns=CHEMBL_ASSAY_DESCRIPTION_REMAPPING)
    pivoted = pivoted.groupby(level=0, axis=1).mean()

    # Normalise data, set "nearly 0" to zero for convenience
    pivoted = (pivoted - pivoted.min()) / (pivoted.max() - pivoted.min())

    # Only select molecules with >2 measurements each (or >n now)
    pivoted = pivoted[pivoted.notna().sum(axis=1).gt(min_measurements_per_mol)]
    # Fill nans with mean
    pivoted = pivoted.fillna(value=pivoted.mean())

    # subtract mean to further normalize then set almost 0 values to 0
    pivoted = pivoted - pivoted.mean()
    pivoted[np.abs(pivoted) < atol] = 0

    # Some sets of assays (e.g. HPLC) are orthogonal to others, so drop these columns
    pivoted = pivoted.dropna(axis=1)

    if verbose:
        print(f"Began with {len(chembl_a_p_raw)} rows, for {len(set(chembl_a_p_raw.description.to_list()))} assays and {len(set(chembl_a_p_raw.chembl_id.to_list()))} molecules.\n")
        print(f"Filtered to {len(chembl_a_p)} rows, for {len(set(chembl_a_p.description.to_list()))} assays and {len(set(chembl_a_p.chembl_id.to_list()))} moleucles.\n")
        print(f"{chembl_a_p["description"].value_counts()}")


    return pivoted
