import sys
import numpy as np
from astropy.table import Table, Column

# Don't print Astropy warnings
import warnings
from astropy.utils.exceptions import AstropyWarning

warnings.filterwarnings("ignore", category=AstropyWarning)


def recioblanco22_calibration(
    uncaldata,
    funcdata,
    val,
    table=None,
    extrapolflag=0,
    func="logg",
    flagcol="FLAGS_GSPSPEC_GAIA",
    floatfill=-9999,
):
    """
    Table 3 and 4 of Recio-Blanco+202 - Analysis of RVS spectra by the General Stellar Parametrizer from spectroscopy.
    Important note! Note that the trends in metallicity and so on are functions of *uncalibrated* log(g).

    uncaldata : numpy array
        Array containing the uncalibrated values to be calibrated.
    funcdata : numpy array
        Array containing the dimension the calibration is parameterised over, e.g. log(g) or Teff
    val : str
        Name of the column to calibrate
    table : astropy table or None
        If not None, the requirement on `extrapol` found in the tables will be folded in with the recommended range in `validflag`.
    extrapolflag : int
        For [alpha/Fe] and [Ca/Fe], two different functions as a function of log(g) is provided depending on the extrapol flag being either 0 og less than or equal
    func : str
        can be either 'logg' or 'teff' and determines whether to use the calibration as a function of either log(g) or t=Teff/5750 K.

    Returns
    -------
    caldata : numpy array
        The calibrated data
    validflag : numpy array
        A bool flag of whether the calibrated values are within the recommended interval or not.
    """
    assert func in ["teff", "logg"]
    if table is not None:
        assert flagcol in table.columns

    if val in ["LOGG"]:
        assert func == "logg"
        ps = [0.4496, -0.0036, -0.0224, 0, 0]
        fmin = 0
        fmax = 10
        maxextrapol = 9
    elif val in ["MH", "MEH"]:
        assert func == "logg"
        ps = [0.274, -0.1373, -0.0050, 0.0048, 0]
        fmin = 0
        fmax = 10
        maxextrapol = 9
    elif val in ["MH_OC"]:
        assert func == "logg"
        ps = [-0.7541, 1.8108, -1.1779, 0.2809, -0.0222]
        fmin = 0
        fmax = 10
        maxextrapol = 9
    elif val in ["ALPHAFE"]:
        if func == "logg":
            if extrapolflag == 0:
                ps = [-0.5809, 0.7018, -0.2402, 0.0239, 0.0000]
                fmin = 1.01
                fmax = 4.85
                maxextrapol = 0
            else:
                ps = [-0.2838, 0.3713, -0.1236, 0.0106, 0.0002]
                fmin = 0.84
                fmax = 4.44
                maxextrapol = 1
        elif func == "teff":
            ps = [-6.6960, 20.8770, -21.0976, 6.8313, 0.0000]
            fmin = 4000
            fmax = 6830
            maxextrapol = 1
    elif val in ["CAFE"]:
        if func == "logg":
            if extrapolflag == 0:
                ps = [-0.6250, 0.7558, -0.2581, 0.0256, 0.0000]
                fmin = 1.01
                fmax = 4.85
                maxextrapol = 0
            else:
                ps = [-0.3128, 0.3587, -0.0816, -0.0066, 0.0020]
                fmin = 0.84
                fmax = 4.98
                maxextrapol = 1
        elif func == "teff":
            ps = [-7.4577, 23.2759, -23.6621, 7.7657, 0.0000]
            fmin = 4000
            fmax = 6830
            maxextrapol = 1
    elif val in ["MGFE"]:
        assert func == "logg"
        ps = [-0.7244, 0.3779, -0.0421, -0.0038, 0.0000]
        fmin = 1.30
        fmax = 4.38
        maxextrapol = 0
    elif val in ["SFE"]:
        if func == "logg":
            ps = [-17.6080, 12.3239, -2.8595, 0.2192, 0.0000]
            fmin = 3.38
            fmax = 4.81
            maxextrapol = 0
        elif func == "teff":
            ps = [0.1930, -0.2234, 0.0000, 0.0000, 0.0000]
            fmin = 5700
            fmax = 6800
            maxextrapol = 0
    elif val in ["SIFE"]:
        assert func == "logg"
        ps = [-0.3491, 0.3757, -0.1051, 0.0092, 0.0000]
        fmin = 1.28
        fmax = 4.85
        maxextrapol = 0
    elif val in ["TIFE"]:
        assert func == "logg"
        ps = [-0.2656, 0.4551, -0.1901, 0.0209, 0.0000]
        fmin = 1.01
        fmax = 4.39
        maxextrapol = 0
    elif val in ["CRFE"]:
        assert func == "logg"
        ps = [-0.0769, -0.1299, 0.1009, -0.0200, 0.0000]
        fmin = 1.01
        fmax = 4.45
        maxextrapol = 0
    elif val in ["FEIH"]:
        assert func == "logg"
        ps = [0.3699, -0.0680, 0.0028, -0.0004, 0.0000]
        fmin = 1.01
        fmax = 4.85
        maxextrapol = 0
    elif val in ["FEIIH"]:
        assert func == "logg"
        ps = [35.5994, -27.9179, 7.1822, -0.6086, 0.0000]
        fmin = 3.53
        fmax = 4.82
        maxextrapol = 0
    elif val in ["NIFE"]:
        assert func == "logg"
        ps = [-0.2902, 0.4066, -0.1313, 0.0105, 0.0000]
        fmin = 1.41
        fmax = 4.81
        maxextrapol = 0
    elif val in ["NFE"]:
        assert func == "logg"
        ps = [0.0975, -0.0293, 0.0238, -0.0071, 0.0000]
        fmin = 1.21
        fmax = 4.79
        maxextrapol = 0

    if func == "teff":
        b = funcdata / 5750
    else:
        b = funcdata

    floatmask = uncaldata == floatfill
    caldata = uncaldata + sum([p * b**i for i, p in enumerate(ps)])
    caldata[floatmask] = floatfill

    validflag = np.ones(len(caldata), dtype=bool)
    validflag[(funcdata < fmin) | (funcdata > fmax)] = False

    if table is not None:
        m = flags_gspspec_check(
            generate_flagpositions(["extrapol"]),
            [maxextrapol],
            operations=["le"],
            table=table,
            col=flagcol,
        )
        validflag |= m

    return caldata, validflag


def calibrate_table(table: Table, col: str, func: str = "logg"):
    if func == "teff":
        funccol = "TEFF_GSPSPEC_GAIA"
        funcstr = "_TEFF"
    else:
        funccol = "LOGG_GSPSPEC_GAIA"
        funcstr = "_LOGG"

    uncalstr = "_UNCALIBRATED"
    calstr = "_CALIBRATED"
    validcol = "FLAG_CAL_" + col + funcstr
    newcolname = col + calstr + funcstr
    uncalcolname = col + uncalstr

    gspspeccol2key = {
        "ALPHAFE_GSPSPEC_GAIA": "ALPHAFE",
        "CAFE_GSPSPEC_GAIA": "CAFE",
        "CRFE_GSPSPEC_GAIA": "CRFE",
        "LOGG_GSPSPEC_GAIA": "LOGG",
        "MGFE_GSPSPEC_GAIA": "MGFE",
        "MH_GSPSPEC_GAIA": "MH",
        "NFE_GSPSPEC_GAIA": "NFE",
        "NIFE_GSPSPEC_GAIA": "NIFE",
        "SFE_GSPSPEC_GAIA": "SFE",
        "SIFE_GSPSPEC_GAIA": "SIFE",
        "TIFE_GSPSPEC_GAIA": "TIFE",
        "FEM_GSPSPEC_GAIA": "FEIH",
        "FEIIM_GSPSPEC_GAIA": "FEIIH",
    }

    val = gspspeccol2key[col]

    if uncalcolname in table.columns:
        uncaldata = table[uncalcolname]
    elif col in table.columns:
        table.rename_column(col, uncalcolname)
        uncaldata = table[uncalcolname]
    else:
        raise ValueError(f"{col}")

    if funccol in table.columns:
        funcdata = table[funccol]
    elif funccol + uncalstr in table.columns:
        funcdata = table[funccol + uncalstr]
    else:
        raise ValueError(f"{funccol}")

    caldata, flagdata = recioblanco22_calibration(
        uncaldata,
        funcdata,
        val=val,
        table=table,
        extrapolflag=0,
        func=func,
        flagcol="FLAGS_GSPSPEC_GAIA",
    )

    # Add columns
    if newcolname in table.columns:
        table.remove_column(newcolname)
    if validcol in table.columns:
        table.remove_column(validcol)

    table.add_column(caldata, name=newcolname)
    table.add_column(flagdata, name=validcol)

    return table


def calibrate_all_gspspec(table: Table, verbose: bool = False):
    calcols = [
        "ALPHAFE_GSPSPEC_GAIA",
        "CAFE_GSPSPEC_GAIA",
        "CRFE_GSPSPEC_GAIA",
        "LOGG_GSPSPEC_GAIA",
        "MGFE_GSPSPEC_GAIA",
        "MH_GSPSPEC_GAIA",
        "NFE_GSPSPEC_GAIA",
        "NIFE_GSPSPEC_GAIA",
        "SFE_GSPSPEC_GAIA",
        "SIFE_GSPSPEC_GAIA",
        "TIFE_GSPSPEC_GAIA",
        "FEM_GSPSPEC_GAIA",
        "FEIIM_GSPSPEC_GAIA",
    ]

    func = "logg"
    for col in calcols:
        if col in table.columns:
            table = calibrate_table(table, col, func=func)
        else:
            if verbose:
                print(f"{col} not in table for {func} calibration")

    calcols = [
        "ALPHAFE_GSPSPEC_GAIA",
        "CAFE_GSPSPEC_GAIA",
        "SFE_GSPSPEC_GAIA",
    ]

    func = "teff"
    for col in calcols:
        if col in table.columns:
            table = calibrate_table(table, col, func=func)
        else:
            print(f"{col} not in table for {func} calibration")
    return table


def add_feh(table, cal=True, floatfill=-9999):
    if cal:
        mh = "MH_GSPSPEC_GAIA_CALIBRATED_LOGG"
        fem = "FEM_GSPSPEC_GAIA_CALIBRATED_LOGG"
        fehname = "FEH_GSPSPEC_GAIA_CALIBRATED_LOGG"
    else:
        mh = "MH_GSPSPEC_GAIA_UNCALIBRATED"
        fem = "FEM_GSPSPEC_GAIA_UNCALIBRATED"
        fehname = "FEH_GSPSPEC_GAIA_UNCALIBRATED"
    if fehname in table.columns:
        table.remove_column(fehname)
    m = (table[mh] == floatfill) | (table[fem] == floatfill)
    feh = table[fem] + table[mh]
    feh[m] = floatfill

    table.add_column(feh, name=fehname)
    return table
