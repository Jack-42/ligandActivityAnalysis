"""
@author Jack Ringer
Date: 6/16/2025
Description:
From a given list of assays, get all compounds which
meet the given criteria for being "active".
"""

from psycopg2 import sql
from psycopg2.extras import DictCursor


def get_active_mols(
    cursor: DictCursor, assay_ids: list[int], pchembl_min_value: int
) -> list[int]:
    query = sql.SQL(
        """
        SELECT md.molregno,md.chembl_id,a.assay_id,a.type,a.value,a.units,a.pchembl_value
        FROM activities a
        JOIN molecule_dictionary md on a.molregno = md.molregno
        WHERE a.assay_id in {assay_ids}
        AND a.pchembl_value >= {pchembl_min_value}
        AND md.structure_type = 'MOL'
        AND (a.data_validity_comment is NULL OR a.data_validity_comment = 'Manually validated';
        """
    ).format(
        assay_ids=sql.Literal(tuple(assay_ids)),
        pchembl_min_value=sql.Literal(pchembl_min_value),
    )
    cursor.execute(query)
    result = cursor.fetchall()
    return result


def main():
    pass


if __name__ == "__main__":
    main()
