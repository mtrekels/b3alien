# b3alien: a Python package to calculate the Target 6.1 headline indicator of the CBD

[![Build](https://github.com/mtrekels/b3alien/actions/workflows/test.yml/badge.svg)](https://github.com/mtrekels/b3alien/actions/workflows/test.yml)
[![Coverage](https://codecov.io/gh/mtrekels/b3alien/branch/main/graph/badge.svg)](https://codecov.io/gh/mtrekels/b3alien)
[![PyPI](https://img.shields.io/pypi/v/b3alien.svg)](https://pypi.org/project/b3alien/)
[![Docs](https://readthedocs.org/projects/b3alien/badge/?version=latest)](https://b3alien.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/mtrekels/b3alien.svg)](https://github.com/mtrekels/b3alien/blob/main/LICENSE)

## Introduction

The historic Kunming-Montreal Global Biodiversity
Framework, which supports the achievement of the
Sustainable Development Goals and builds on the
Convention on Biological Diversity’s (CBD) previous
Strategic Plans, sets out an ambitious pathway to reach
the global vision of a world living in harmony with nature
by 2050. Among the Framework’s key elements are 23
targets for 2030. In order to track the progress on the
targets, a number of indicators were agreed upon for
each target. The B3ALIEN software provides a technical
solution to track Target 6: “Reduce the Introduction of
Invasive Alien Species by 50% and Minimize Their
Impact.” It mainly focusses on the headline indicator: rate
of invasive alien species establishment, but can provide
input to some of the complementary indicators.

Decision makers at local, regional, national and
international levels need accurate and reliable
information about status, trends, threats, and they need
data presented in an actionable and understandable
format, with measures of uncertainty. Furthermore, we
need synthesized data products that can be combined
with other environmental data, such as climate, soil
chemistry, land use, altitude... B3ALIEN is built upon the
concept of data cubes developed in the Horizon Europe
Biodiversity Building Blocks for Policy project (b-
cubed.eu). It uses the solid foundations of the GBIF
infrastructure, where tools such as the GBIF Taxonomic
Backbone and the Global Registry of Introduced and
Invasive Species are available by default. Readily available occurrence data is used to determine and estimate
accurately the rate of introduction of alien species..

## Architecture

![png](docs/_static/images/architecture_b3alien.png)

## Example usage

```python
from b3alien import b3cube
from b3alien import visualisation
from b3alien.utils.runtime import in_jupyter
from b3alien import griis
from b3alien import simulation

import matplotlib
matplotlib.use("TkAgg")
```



```python
cube = b3cube.OccurrenceCube("gs://b-cubed-eu/data_PT-30b.parquet", gproject='$GPROJECT-ID')
```



```python
print(cube.df)
```



           kingdom  kingdomkey           phylum  phylumkey              class  \
    0      Plantae           6       Charophyta    7819616       Charophyceae   
    1      Plantae           6     Tracheophyta    7707728     Polypodiopsida   
    2      Plantae           6     Tracheophyta    7707728     Polypodiopsida   
    3      Plantae           6     Tracheophyta    7707728     Polypodiopsida   
    4      Plantae           6     Tracheophyta    7707728     Polypodiopsida   
    ...        ...         ...              ...        ...                ...   
    52187  Plantae           6     Tracheophyta    7707728      Magnoliopsida   
    52188  Plantae           6     Tracheophyta    7707728      Magnoliopsida   
    52189  Plantae           6     Tracheophyta    7707728      Magnoliopsida   
    52190  Plantae           6  Marchantiophyta          9  Jungermanniopsida   
    52191  Plantae           6        Bryophyta         35          Bryopsida   
    
           classkey            order  orderkey          family  familykey  ...  \
    0           328         Charales       626       Characeae       8782  ...   
    1       7228684      Salviniales   7229405    Salviniaceae       6629  ...   
    2       7228684     Polypodiales       392   Polypodiaceae       2368  ...   
    3       7228684     Polypodiales       392   Polypodiaceae       2368  ...   
    4       7228684     Polypodiales       392   Polypodiaceae       2368  ...   
    ...         ...              ...       ...             ...        ...  ...   
    52187       220     Malpighiales      1414   Euphorbiaceae       4691  ...   
    52188       220     Malpighiales      1414   Euphorbiaceae       4691  ...   
    52189       220     Malpighiales      1414   Euphorbiaceae       4691  ...   
    52190       126  Jungermanniales       381  Calypogeiaceae       2289  ...   
    52191       327        Pottiales       621      Pottiaceae       4671  ...   
    
          classcount  ordercount familycount genuscount distinctobservers  \
    0              1           1           1          1                 1   
    1             11           1           1          1                 1   
    2              1           1           1          1                 1   
    3              1           1           1          1                 1   
    4              6           6           1          1                 1   
    ...          ...         ...         ...        ...               ...   
    52187         27           4           4          2                 1   
    52188         27           4           4          2                 1   
    52189         26           3           3          3                 1   
    52190          2           2           1          1                 1   
    52191          4           1           1          1                 1   
    
          occurrences  mintemporaluncertainty  mincoordinateuncertaintyinmeters  \
    0               1                 2678400                            1000.0   
    1               1                      60                              23.0   
    2               1                   86400                            1000.0   
    3               1                 2505600                            1000.0   
    4               1                   86400                            1000.0   
    ...           ...                     ...                               ...   
    52187           1                       1                              16.0   
    52188           1                       1                               6.0   
    52189           3                   86400                            2400.0   
    52190           1                       1                            1000.0   
    52191           1                   86400                              30.0   
    
               cellCode                                           geometry  
    0      W016N32ACAAA  POLYGON ((-17 32.71875, -16.96875 32.71875, -1...  
    1      W016N32ACADB  POLYGON ((-16.90625 32.65625, -16.875 32.65625...  
    2      W016N32AACCA  POLYGON ((-17 32.78125, -16.96875 32.78125, -1...  
    3      W016N32AACCD  POLYGON ((-16.96875 32.75, -16.9375 32.75, -16...  
    4      W016N32AACCD  POLYGON ((-16.96875 32.75, -16.9375 32.75, -16...  
    ...             ...                                                ...  
    52187  W017N32BDBCB  POLYGON ((-17.09375 32.65625, -17.0625 32.6562...  
    52188  W017N32BDBCB  POLYGON ((-17.09375 32.65625, -17.0625 32.6562...  
    52189  W017N32BDBDA  POLYGON ((-17.0625 32.65625, -17.03125 32.6562...  
    52190  W016N32AACDB  POLYGON ((-16.90625 32.78125, -16.875 32.78125...  
    52191  W017N32BBCAD  POLYGON ((-17.21875 32.8125, -17.1875 32.8125,...  
    
    [52192 rows x 28 columns]



```python
cube._species_richness()
```





```python
print(cube.richness)
```




                 cell  richness
    0    W016N30DDDDA         3
    1    W016N32AACAC       212
    2    W016N32AACAD       152
    3    W016N32AACBC       168
    4    W016N32AACBD       225
    ..            ...       ...
    120  W017N32BDBCA       136
    121  W017N32BDBCB       165
    122  W017N32BDBDA       189
    123  W017N32BDBDB       276
    124  W017N32BDBDD        91
    
    [125 rows x 2 columns]



```python
b3cube.plot_richness(cube.richness, cube.df)
```

![png](docs/_static/images/richness_plot.png)


```python
CL = griis.CheckList("$YOUR_DIRECTORY/merged_distr.txt")
```


```python
d_s, d_c = b3cube.cumulative_species(cube, CL.species)
```

```python
time, rate = b3cube.calculate_rate(d_c)
```


```python
C1 = simulation.simulate_solow_costello(time, rate, vis=True)
```



    Optimization terminated successfully.
             Current function value: -263.092115
             Iterations: 172
             Function evaluations: 287


    
![png](docs/_static/images/output_9_2.png)
    


```python

```