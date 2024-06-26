# This is a brief collection of functions from Skycats.
import io
import os
import time
import sys
import ssl
import csv
import traceback
import warnings
import numpy as np
from natsort import natsorted
from astropy.table import Table, join, vstack, unique, Column
from astroquery.utils.tap.core import TapPlus

from zero_point import zpt


# Gaia specific
def find_gaia_from_sourceids(
    filename: str,
    keytable: Table,
    dr: str = "DR3",
    gaiacat: str | None = None,
    source_id: str = "source_id",
    frmt: str = "votable",
    cols: list | None = None,
    chunk: int = 1000,
):
    if os.path.exists(filename):
        table = Table.read(filename, format=frmt)
        table = standardize(table)
    else:
        chunks = []
        gaia = TapPlus(url="https://gea.esac.esa.int/tap-server/tap")

        if dr == "DR2":
            if not "SOURCE_ID_GAIADR2" in keytable.columns:
                keytable.rename_column("SOURCE_ID_GAIA", "SOURCE_ID_GAIADR2")
            stars = keytable["SOURCE_ID_GAIADR2"].data.astype(str)
            if gaiacat is None:
                gaiacat = "gaiadr2.gaia_source"
        elif dr == "eDR3":
            if not "SOURCE_ID_GAIAEDR3" in keytable.columns:
                keytable.rename_column("SOURCE_ID_GAIA", "SOURCE_ID_GAIAEDR3")
            stars = keytable["SOURCE_ID_GAIAEDR3"].data.astype(str)
            if gaiacat is None:
                gaiacat = "gaiaedr3.gaia_source"
        elif dr == "DR3":
            stars = keytable["SOURCE_ID_GAIA"].data.astype(str)
            if gaiacat is None:
                gaiacat = "gaiadr3.gaia_source"
        else:
            raise NotImplementedError(dr)

        if cols is None:
            cols = "*"

        # Fix bug in astropy
        # tapContext ends with a slash, but it shouldn't!
        gaia._Tap__connHandler._TapConn__tapContext = "/tap-server/tap"
        # If you get a SSL: DL KEY TOO SMALL error on your server/laptop do this:
        ssl._DEFAULT_CIPHERS += ":!DH"

        for j in range(0, len(stars), chunk):
            print(f"Make {filename}, progress: {j}/{len(stars)}")
            job = gaia.launch_job(
                f"""
                                  SELECT {cols}
                                  FROM {gaiacat} AS cat
                                  WHERE cat."{source_id}"
                IN ({",".join(adql_escape(s) for s in stars[j : j + chunk])})
                """,
                verbose=False,
            )
            tempttable = job.get_results()
            chunks.append(tempttable)
        table = vstack(chunks)

        if "gaia_source" in gaiacat:
            table = standardize(table, catid="GAIA")
            if dr == "DR2":
                table.rename_column("SOURCE_ID_GAIA", "SOURCE_ID_GAIADR2")
            elif dr == "eDR3":
                table.rename_column("SOURCE_ID_GAIA", "SOURCE_ID_GAIAEDR3")
        elif "astrophysical_parameters" in gaiacat:
            table = standardize(table, catid="GAIA", debug=True)
        elif "dr2_neighbourhood" in gaiacat:
            table = standardize(table)
        else:
            table = standardize(table, catid="BAILERJONES")

        table.write(filename, format=frmt)
    return table


def compute_gaiaphotometryerror(
    table: Table,
    dr: str = "DR3",
    verbose: bool = False,
    error_gmag: str = "ERROR_PHOT_G_MEAN_FLUX_GAIA",
    gmag: str = "phot_g_mean_flux_GAIA".upper(),
    error_bpmag: str = "ERROR_PHOT_BP_MEAN_FLUX_GAIA",
    bpmag: str = "phot_bp_mean_flux_GAIA".upper(),
    error_rpmag: str = "ERROR_PHOT_RP_MEAN_FLUX_GAIA",
    rpmag: str = "phot_rp_mean_flux_GAIA".upper(),
):
    # Compute photometry errors in magnitudes and a floor of 2.3 mmag in
    # quadrature as discussed in Arenou et al. 2018
    usedcols = [gmag, error_gmag]
    for uc in usedcols:
        if uc not in table.columns:
            if verbose:
                print(f"The column {uc} is not in the table")
                print("SKYCATS skips computing photometric uncertainties in Gaia\n")
            return table

    def magerror(error_flux, flux, table, relcoeff, sys):
        assert error_flux in table.columns
        assert flux in table.columns
        error_mag = np.sqrt(
            (relcoeff * (table[error_flux].data / table[flux].data)) ** 2 + sys**2
        )
        return error_mag

    if dr in ["eDR3", "DR3"]:
        # Following Riello+2021, the are no systematic trends larger than 1 mmag mag-1.
        # See their Fig 25.
        # https://vizier.cds.unistra.fr/viz-bin/VizieR-n?-source=METAnot&catid=1350&notid=63&-out=text
        relcoeff = 1
        G_err = magerror(
            error_gmag,
            gmag,
            table,
            relcoeff,
            0.0027553202,
        )
        BP_err = magerror(
            error_bpmag,
            bpmag,
            table,
            relcoeff,
            0.0027901700,
        )
        RP_err = magerror(
            error_rpmag,
            rpmag,
            table,
            relcoeff,
            0.0037793818,
        )
    else:
        print("Gaia data release version unknown\n")
        return table

    G_err_col = Column(G_err, name="ERROR_PHOT_G_MEAN_MAG_GAIA")
    BP_err_col = Column(BP_err, name="ERROR_PHOT_BP_MEAN_MAG_GAIA")
    RP_err_col = Column(RP_err, name="ERROR_PHOT_RP_MEAN_MAG_GAIA")
    table.add_column(G_err_col)
    table.add_column(BP_err_col)
    table.add_column(RP_err_col)
    return table


def adql_escape(star: str):
    """
    A function that can make stings in the right format for adql queries.
    Inspiration from https://github.com/
    gmantele/taplib/blob/master/src/adql/query/operand/StringConstant.java
    """
    return "'" + star.replace("'", "''") + "'"


def add_parallaxwithoffset(
    table: Table,
    dr: str = "DR3",
    parallaxcol: str = "parallax_GAIA".upper(),
    singleoffset: bool = False,
    verbose: bool = False,
    gmag: str = "phot_g_mean_mag_GAIA".upper(),
    pseudocolour: str = "pseudocolour_GAIA".upper(),
    nueff: str = "NU_EFF_ASTROMETRY_GAIA",
    ecl: str = "ecl_lat_GAIA".upper(),
    paramssolved: str = "astrometric_params_solved_GAIA".upper(),
):
    # Here we add the parallax columns with different offsets:
    assert parallaxcol in table.columns, list(table.columns)
    uncorrparallaxcolname = f"UNCORRECTED_{parallaxcol.upper()}"
    table.rename_column(parallaxcol, uncorrparallaxcolname)

    usedcols = [uncorrparallaxcolname, gmag, pseudocolour, nueff, ecl, paramssolved]
    for uc in usedcols:
        if uc not in table.columns:
            if verbose:
                print(f"The column {uc} is not in the table")
                print("SKYCATS skips adding offsets to Gaia parallax\n")
            return table

    if dr in ["eDR3", "DR3"]:
        if singleoffset:
            # Gaia Collaboration 2021a (https://arxiv.org/abs/2012.01533) offset
            # of 17 mas computed from quasars.
            floatfill = get_fillvalue(table[uncorrparallaxcolname])
            corrpar = np.ones(len(table)) * floatfill
            corrparmask = table[uncorrparallaxcolname] != floatfill

            st_par = np.copy(corrpar)
            st_par[corrparmask] = (
                table[uncorrparallaxcolname][corrparmask] + 0.017
            )
            st_par = Column(st_par, name="PARALLAX_GLOBALOFFSET_GAIA")
            table.add_column(st_par)
        else:
            # Here I will use the correction from
            # https://gitlab.com/icc-ub/public/gaiadr3_zeropoint
            usedcols = [gmag, pseudocolour, nueff, ecl, paramssolved]
            for uc in usedcols:
                if uc not in table.columns:
                    if verbose:
                        print(f"The column {uc} is not in the table")
                        print("SKYCATS skips adding offsets to Gaia parallax\n")
                    return table

            zpt.load_tables()

            floatfill = get_fillvalue(table[uncorrparallaxcolname])
            corrpar = np.ones(len(table)) * floatfill
            corrparmask = table[uncorrparallaxcolname] != floatfill

            corr = zpt.get_zpt(
                table[gmag][corrparmask],
                table[nueff][corrparmask],
                table[pseudocolour][corrparmask],
                table[ecl][corrparmask],
                table[paramssolved][corrparmask],
                _warnings=True,
            )
            corr[~np.isfinite(corr)] = 0.0

            # As explained in Lindegren et al. 2020, the interpolations are only
            # calibrated within the following intervals - outside it is extrapolation
            warnings = (
                (6 < table[gmag][corrparmask]) & (table[gmag][corrparmask] < 21) &
                # 5-param solutions
                (
                    (table[paramssolved][corrparmask] == 31)
                    & (1.1 < table[nueff][corrparmask])
                    & (table[nueff][corrparmask] < 1.9)
                )
                |
                # 6-param solutions
                (
                    (table[paramssolved][corrparmask] == 95)
                    & (1.24 < table[pseudocolour][corrparmask])
                    & (table[pseudocolour][corrparmask] < 1.72)
                )
            )

            # Add warnings column
            warncol = np.zeros(len(table), dtype=bool)
            warncol[corrparmask] = ~warnings
            warncol = Column(warncol, name="WARNINGS_PARALLAXOFFSET_GAIA")

            # Add corrected parallax as a column
            st_par = np.copy(corrpar)
            st_par[corrparmask] = (
                table[uncorrparallaxcolname][corrparmask] - corr
            )
            st_par = Column(st_par, name="PARALLAX_GAIA")
            table.add_column(st_par)
            table.add_column(warncol)
    else:
        print("Unknown Gaia data release version")

    return table


# Clean and groom cats
def standardize(
    table: Table,
    catid: str | None = None,
    label: None | str = None,
    old: str | None = None,
    debug: bool = False,
    remove: bool = True,
):
    # This function is a very scrapped version of the one in Skycats.
    table = check_columns_for_objects(table)
    for col in table.columns:
        fillvalue = get_fillvalue(table[col])
        table[col].fill_value = fillvalue
        table[col].description = ""
        table[col].unit = None
    table = table.filled()
    if catid is not None:
        table = rename_columns(table, catid, debug=debug)
    table = make_ids_to_string(table)
    table = make_objects_to_unicode(table)
    table = replace_old_fillvalue(table, old)
    return table


the_rename_mapping = None


def get_rename_mapping():
    #global the_rename_mapping
    #if the_rename_mapping is not None:
    #    return the_rename_mapping
    mapping = {}
    here = os.path.realpath(os.path.dirname(__file__))
    with open(os.path.join(here, "column_overview.txt"), newline="") as fp:
        read = csv.reader(fp)
        for line in read:
            column = line[0]
            oldnames = line[1].split()
            for oldname in oldnames:
                mapping[oldname] = column
    the_rename_mapping = mapping
    return mapping


def rename_columns(table, catid, debug=False):
    colnames = table.colnames
    mapping = get_rename_mapping()
    if debug:
        for col in colnames:
            print(
                f"{col.upper()}_{catid.upper()},{col}_{catid},category,definition from {catid},0,0"
            )
    for col in colnames:
        if col.endswith(str(catid)):
            COL = col.replace(" ", "")
        else:
            COL = col.replace(" ", "") + "_" + str(catid)
        if not COL in list(mapping.values()):
            new_name = mapping[COL]
            table.rename_column(col, new_name)
    return table


def check_columns_for_objects(table: Table, verbose: bool = False):
    for col in table.colnames:
        if table[col].dtype == "object":
            if verbose:
                print(f"Object found! Column {str(col)}")
            table[col] = table[col].astype(str)
    return table


def get_fillvalue(col):
    if col.dtype.kind in ["S", "U"]:
        return ""
    elif col.dtype.kind in ["i", "f", "b", "u"]:
        return -9999.0
    else:
        raise NotImplementedError(col.dtype.kind)


def make_ids_to_string(table: Table):
    ids = [
        "LABEL",
        "KIC_ID",
        "TMASS_ID",
        "HIP_ID",
        "EPIC_ID",
        "HIP_ID",
        "TYC_ID",
        "UCAC4_ID",
        "SDSS_ID",
        "APOGEE_ID",
        "RAVE_ID",
        "KOI_ID",
        "KIS_ID",
        "SOURCE_ID_GAIA",
        "SOURCE_ID_GAIADR2",
        "SOURCE_ID_GAIAEDR3",
        "TIC_ID",
    ]
    for id in ids:
        if id in table.columns:
            table[id] = table[id].astype("U")
    return table


def make_objects_to_unicode(table: Table):
    for col in table.columns:
        if table[col].dtype.kind in ["S", "object"]:
            table[col] = table[col].astype("str")
    return table


def replace_old_fillvalue(table: Table, olds: list | float | str):
    if type(olds) is not list:
        olds = [olds]
    for old in olds:
        for col in table.columns:
            new = get_fillvalue(table[col])
            if old is None:
                old = new
            if isinstance(old, str):
                mask = table[col].data.astype("U") == old
            else:
                if np.isnan(old):
                    if table[col].dtype.kind not in ["S", "U", "object"]:
                        mask = np.isnan(table[col].data)
                    else:
                        mask = table[col].data == "nan"
                else:
                    if table[col].dtype.kind not in ["S", "U", "object"]:
                        mask = table[col].data == old
                    else:
                        mask = np.zeros(len(table[col].data), dtype=bool)
            table[col][mask] = new
    return table


def delete_unusedcols(table):
    cols = []
    ignore = []
    co = os.path.join(os.path.dirname(__file__), "column_overview.txt")
    with open(co, newline="") as fp:
        read = csv.reader(fp)
        for line in read:
            colname = line[0]
            category = line[2]
            notused = line[4]
            if category in ["DELETE"]:
                ignore.append(colname)
            elif notused == "1":
                ignore.append(colname)
            else:
                cols.append(colname)

    cols = set(cols)
    ignore = set(ignore)
    tablecols = set(table.columns)
    unknowncols = tablecols - cols - ignore
    if unknowncols:
        print("Unknown columns:")
        for c in sorted(unknowncols):
            print(f"{c},")
        raise SystemExit

    catcolumns = sorted(cols & tablecols)
    table = table[catcolumns]
    return table
