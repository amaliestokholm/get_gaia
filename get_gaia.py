# Get Gaia values from Gaia DR3 source_id
# Amalie Stokholm

import numpy as np
from astropy.table import Table, join
import funkykitten
import gspspec

# Specify stuff here
starid = "SOURCE_ID_GAIA"

# Which columns in 'gaiadr3.astrophysical_parameters' do you want to keep?
gspspeccols = [
    "source_id",
    "teff_gspspec",
    "teff_gspspec_lower",
    "teff_gspspec_upper",
    "mh_gspspec",
    "mh_gspspec_lower",
    "mh_gspspec_upper",
    #    "sife_gspspec",
    #    "sife_gspspec_lower",
    #    "sife_gspspec_upper",
    "alphafe_gspspec",
    "alphafe_gspspec_lower",
    "alphafe_gspspec_upper",
    "logg_gspspec",
    "logg_gspspec_lower",
    "logg_gspspec_upper",
    "flags_gspspec",
]


# Actually do stuff
a = Table.read(targetlist, format="ascii.commented_header")
print(f"Length of table is {len(a)}")
print(f"Number of unique targets is {len(np.unique(a[starid]))}")
if starid == "sourceid":
    a.rename_column("sourceid", starid)

if not os.path.exists(tempdir):
    os.makedirs(tempdir)

if isinstance(a, Column):
    a = Table(data=a, name=starid)

gaiatable = funkykitten.find_gaia_from_sourceids(
    os.path.join(tempdir, f"{projectname}_gaia.vot"),
    keytable=a,
    dr="DR3",
)
gaiatable = funkykitten.standardize(gaiatable)

if "source_id" in list(gaiatable.columns):
    gaiatable.rename_column("source_id", starid)

gaiatable = funkykitten.add_parallaxwithoffset(gaiatable, dr="DR3")
gaiatable.rename_column("PARALLAX_GAIA", "ZPCORR_PARALLAX_GAIA")

gspspeccols = ", ".join(gspspeccols)

gspspectable = funkykitten.find_gaia_from_sourceids(
    os.path.join(tempdir, f"{projectname}_gspspecgaia.vot"),
    keytable=a,
    dr="DR3",
    gaiacat="gaiadr3.astrophysical_parameters",
    cols=gspspeccols,
)

gspspectable = gspspec.calibrate_all_gspspec(gspspectable)

# Get overview
for col in list(gspspectable.columns):
    if "_CALIBRATED" in col:
        print(col)
    elif col == "source_id":
        gspspectable.rename_column(col, starid)

gaiatable = join(gaiatable, gspspectable, join_type="left", keys=starid)

gaiatable = funkykitten.compute_gaiaphotometryerror(gaiatable, dr="DR3", verbose=True)

gaiatable.write(
    f"./data/gaiatable_{projectname}.ascii",
    format="ascii.commented_header",
    overwrite=True,
)

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

if not os.path.exists(resultdir):
    os.makedirs(resultdir)

resultsfile = os.path.join(resultdir, resultsfile)
print(f"Writing the following columns to {resultsfile}: {tokeep}")

b = join(
    a,
    gaiatable[tokeep],
    join_type="left",
)
b.rename_column(starid, "starid")
b.write(resultsfile, format="ascii.commented_header", overwrite=True)
