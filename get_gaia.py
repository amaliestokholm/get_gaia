# Get Gaia values from Gaia DR3 source_id
# Amalie Stokholm

import os
import argparse
import numpy as np
from astropy.table import Table, Column, join

import funkykitten
import gspspec

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--identifier")
parser.add_argument("-p", "--projectname", default="mygaiaquery")
parser.add_argument("-s", "--starid", default="SOURCE_ID_GAIA")
parser.add_argument("-t", "--targetlist", default="./example/targetlist.ascii")
parser.add_argument("--dr", default="DR3")
parser.add_argument("--tempdir", default="./example/gaia")
parser.add_argument("--resultdir", default="./results")
parser.add_argument("--resultsfile", default="sample.ascii")
parser.add_argument("--gspspec", action="store_true")


def get_gaiamain(
    a,
    projectname: str,
    dr: str = "DR3",
    starid: str = "SOURCE_ID_GAIA",
    tempdir: str = "./example/gaia/",
    resultdir: str = "./results",
    resultsfile: str = "sample.ascii",
):
    gaiatable = funkykitten.find_gaia_from_sourceids(
        os.path.join(tempdir, f"{projectname}_gaia.vot"),
        keytable=a,
        dr=dr,
    )

    if "source_id" in list(gaiatable.columns):
        gaiatable.rename_column("source_id", starid)

    gaiatable = funkykitten.add_parallaxwithoffset(gaiatable, dr=dr)

    gaiatable.rename_column("PARALLAX_GAIA", "ZPCORR_PARALLAX_GAIA")

    gaiatable = funkykitten.compute_gaiaphotometryerror(gaiatable, dr=dr, verbose=True)

    gaiatable.write(
        os.path.join(tempdir, f"gaiatable_{projectname}.ascii"),
        format="ascii.commented_header",
        overwrite=True,
    )

    resultsfile = os.path.join(resultdir, resultsfile)
    print(f"Write results from the Gaia main catalogue to {resultsfile}")
    gaiatable.write(resultsfile, format="ascii.commented_header", overwrite=True)

    return gaiatable


def get_gspspec(
    a: Table,
    gaiatable: Table,
    projectname: str,
    dr: str = "DR3",
    starid: str = "SOURCE_ID_GAIA",
    tempdir: str = "./example/gaia/",
    resultdir: str = "./results",
    resultsfile: str = "sample.ascii",
    verbose: bool = False,
):
    # Which columns in 'gaiadr3.astrophysical_parameters' do you want to keep?
    gspspeccols = [
        "source_id",
        "teff_gspspec",
        "teff_gspspec_lower",
        "teff_gspspec_upper",
        "mh_gspspec",
        "mh_gspspec_lower",
        "mh_gspspec_upper",
        "sife_gspspec",
        "sife_gspspec_lower",
        "sife_gspspec_upper",
        "alphafe_gspspec",
        "alphafe_gspspec_lower",
        "alphafe_gspspec_upper",
        "logg_gspspec",
        "logg_gspspec_lower",
        "logg_gspspec_upper",
        "flags_gspspec",
    ]
    gspspeccols = ", ".join(gspspeccols)

    gspspectable = funkykitten.find_gaia_from_sourceids(
        os.path.join(tempdir, f"{projectname}_gspspecgaia.vot"),
        keytable=a,
        dr=dr,
        gaiacat="gaiadr3.astrophysical_parameters",
        cols=gspspeccols,
    )

    gspspectable = gspspec.calibrate_all_gspspec(gspspectable)

    # Get overview
    if verbose:
        for col in list(gspspectable.columns):
            if "_CALIBRATED" in col:
                print(col)

    if "source_id" in list(gspspectable.columns):
        gspspectable.rename_column(col, starid)

    a[starid] = a[starid].astype(str)
    gaiatable[starid] = gaiatable[starid].astype(str)
    gaiatable = join(gaiatable, gspspectable, join_type="left", keys=starid)

    keepcols = [
        starid,
        "ZPCORR_PARALLAX_GAIA",
        "WARNINGS_PARALLAXOFFSET_GAIA",
        "TEFF_GSPSPEC_GAIA",
        "LOWER_ERROR_TEFF_GSPSPEC_GAIA",
        "UPPER_ERROR_TEFF_GSPSPEC_GAIA",
        "ALPHAFE_GSPSPEC_GAIA_UNCALIBRATED",
        "LOWER_ERROR_ALPHAFE_GSPSPEC_GAIA",
        "UPPER_ERROR_ALPHAFE_GSPSPEC_GAIA",
        "ALPHAFE_GSPSPEC_GAIA_CALIBRATED_LOGG",
        "FLAG_CAL_ALPHAFE_GSPSPEC_GAIA_LOGG",
        "LOGG_GSPSPEC_GAIA_UNCALIBRATED",
        "LOWER_ERROR_LOGG_GSPSPEC_GAIA",
        "UPPER_ERROR_LOGG_GSPSPEC_GAIA",
        "LOGG_GSPSPEC_GAIA_CALIBRATED_LOGG",
        "FLAG_CAL_LOGG_GSPSPEC_GAIA_LOGG",
        "MH_GSPSPEC_GAIA_UNCALIBRATED",
        "UPPER_ERROR_MH_GSPSPEC_GAIA",
        "LOWER_ERROR_MH_GSPSPEC_GAIA",
        "MH_GSPSPEC_GAIA_CALIBRATED_LOGG",
        "FLAG_CAL_MH_GSPSPEC_GAIA_LOGG",
        "FLAGS_GSPSPEC_GAIA",
        "PHOT_G_MEAN_FLUX_GAIA",
        "ERROR_PHOT_G_MEAN_FLUX_GAIA",
        "PHOT_BP_MEAN_FLUX_GAIA",
        "ERROR_PHOT_BP_MEAN_FLUX_GAIA",
        "PHOT_RP_MEAN_FLUX_GAIA",
        "ERROR_PHOT_RP_MEAN_FLUX_GAIA",
    ]

    a[starid] = a[starid].astype(str)
    gaiatable[starid] = gaiatable[starid].astype(str)
    tokeep = list(set(gaiatable.columns) & set(keepcols))

    resultsfile = os.path.join(
        resultdir, f"{resultsfile.split('.')[0]}_gspspec.{resultsfile.split('.')[-1]}"
    )
    print(f"Writing the following columns to {resultsfile}: {tokeep}")

    b = join(
        a,
        gaiatable[tokeep],
        join_type="left",
    )
    b.rename_column(starid, "starid")
    b.write(resultsfile, format="ascii.commented_header", overwrite=True)


def main():
    args = parser.parse_args()

    # Specify stuff here
    starid = args.starid
    projectname = args.projectname
    targetlist = args.targetlist
    tempdir = args.tempdir
    resultdir = args.resultdir
    resultsfile = args.resultsfile
    gspspec = args.gspspec
    dr = args.dr

    # Actually do stuff
    if args.identifier is None:
        a = Table.read(targetlist, format="ascii.commented_header")
    else:
        a = Table(data=[[args.identifier]], names=[args.starid], dtype=[str])

    print(f"Length of table is {len(a)}")
    print(f"Number of unique targets is {len(np.unique(a[starid]))}")

    if not starid in a.columns:
        assert "sourceid" in a.columns
        a.rename_column("sourceid", starid)

    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # Get Gaia Main
    gaiatable = get_gaiamain(
        a, projectname, dr, starid, tempdir, resultdir, resultsfile
    )

    if gspspec:
        if not os.path.exists(resultdir):
            os.makedirs(resultdir)

        # Get GSP-Spec
        gspspectable = get_gspspec(
            a, gaiatable, projectname, dr, starid, tempdir, resultdir, resultsfile
        )


if __name__ == "__main__":
    main()
