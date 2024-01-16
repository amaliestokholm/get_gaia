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
a = Table.read("./data/targetlist.ascii", format="ascii")
print(
    f"Length of table is {len(a)}, number of unique targets is {len(np.unique(a[starid]))}"
)

if starid == "sourceid":
    a.rename_column("starid", "SOURCE_ID_GAIA")

gaiatable = funkykitten.find_gaia_from_sourceids(
    f"./data/gaia/{projectname}_gaia.vot",
    keytable=a,
    dr="DR3",
)
gaiatable = funkykitten.add_parallaxwithoffset(gaiatable, dr="DR3")
gaiatable.rename_column("PARALLAX_GAIA", "ZPCORR_PARALLAX_GAIA")

gspspeccols = ", ".join(gspspeccols)

gspspectable = funkykitten.find_gaia_from_sourceids(
    f"./data/gaia/{projectname}_gspspecgaia.vot",
    keytable=a["SOURCE_ID_GAIA"],
    dr="DR3",
    gaiacat="gaiadr3.astrophysical_parameters",
    cols=gspspeccols,
)

gspspectable = gspspec.calibrate_all_gspspec(gspspectable)
for col in gspspectable.columns:
    if "_CALIBRATED" in col:
        print(col)

gaiatable = join(gaiatable, gspspectable, join_type="left", keys="SOURCE_ID_GAIA")

gaiatable = funkykitten.compute_gaiaphotometryerror(gaiatable, dr="DR3", verbose=True)

gaiatable.write(
    f"./data/gaiatable_{projectname}.ascii",
    format="ascii.commented_header",
    overwrite=True,
)

keepcols = [
    "SOURCE_ID_GAIA",
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

a["SOURCE_ID_GAIA"] = a["SOURCE_ID_GAIA"].astype(str)

b = join(
    a,
    gaiatable[keepcols],
    join_type="left",
)
b.rename_column("SOURCE_ID_GAIA", "starid")
b.write("./data/sample.ascii", format="ascii.commented_header", overwrite=True)
