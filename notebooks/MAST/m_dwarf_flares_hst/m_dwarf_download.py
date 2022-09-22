"""
External script used to download the list of high proper motion M dwarfs.

author: @arjunsavel
"""
from astroquery.utils.tap.core import TapPlus
import pandas as pd


VIZIER_TAP_URL = 'http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap'
viz = TapPlus(url=VIZIER_TAP_URL)

table = "J/AJ/161/63/table3"

job = viz.launch_job_async(
    f"""SELECT TOP 1000 *
    FROM 
        "{table}"
   """,
    output_format="csv",
)
tab = job.get_results()
df = tab.to_pandas()

df = df.sort_values('Jmag').reset_index()


df.to_csv('high_pm_m_dwarfs.csv', index=False)