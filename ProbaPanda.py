import pandas as pd
import matplotlib.pyplot as plt
#url="https://"
url = "https://raw.githubusercontent.com/santono/gpss_f/master/r_12.csv"
df= pd.read_csv(url,
                          index_col=[2,3], 
#                      usecols=[5],
                          parse_dates=True)
#df.head()
#df[(40,1)]
df.dtypes