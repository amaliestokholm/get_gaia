# Make Gaia querying simple
This is just a small collection of scripts that makes Gaia querying easier.
Most of this are snippets from the `skycats` repository.


## I have a single target and I know the Gaia DR3 source\_id.
```
python3 get_gaia.py -i {yoursourceidhere}
```
for instance
```
python3 get_gaia.py -i 2135550755683405952

```
Alternatively, you can add your Gaia DR3 source\_id to a file
```
python3 get_gaia.py --targetlist "./example/targetlist.ascii"
```

## I want the GSP-Spec values for my target!
Add the flag `--gspspec`
```
python3 get_gaia.py --gspspec
```

## I want the FLAME results based on GSP-Spec instead of GSP-Phot!
This I have not added to this repository yet - just write me or make an issue and I'll add it :)


## I have a list of Gaia DR3 source\_ids
Have a look at the `example` directory here.

In the `targetlist.ascii`, I have a file with Gaia DR3 source ids of a single star. 
You edit this so it contains the Gaia DR3 source ids of all your targets.
Now have a look at the `get_gaia.py` script. This is the main file


## How do I easily add a Gaia query to my Python script?
To query Gaia using a list of source ids:
```
from astroquery.utils.tap.core import TapPlus

stars = ['list', 'of', 'gaiasourceids']
cols = ['columns', 'you', 'want']
gaiacat = "gaiadr3.gaia_source"  # or whichever you prefer

# If you are using the gaia module and not TapPlus, ignore this
gaia = TapPlus(url="https://gea.esac.esa.int/tap-server/tap")
gaia._Tap__connHandler._TapConn__tapContext = "/tap-server/tap"

def adql_escape(star):
    """
    A function that can make stings in the right format for adql queries.
    Inspiration from https://github.com/
    gmantele/taplib/blob/master/src/adql/query/operand/StringConstant.java
    """
    return "'" + star.replace("'", "''") + "'"

job = gaia.launch_job(
                f"""
                                  SELECT {cols}
                                  FROM {gaiacat} AS cat
                                  WHERE cat."{source_id}"
                IN ({",".join(adql_escape(s) for s in stars)})
                """,
                verbose=False,
            )
```
