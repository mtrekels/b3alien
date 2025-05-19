Example usage
=============

.. code-block:: python

    from b3alien import b3cube
    from b3alien import visualisation
    from b3alien.utils.runtime import in_jupyter
    from b3alien import griis
    from b3alien import simulation

    import matplotlib
    matplotlib.use("TkAgg")

.. code-block:: python

    cube = b3cube.OccurrenceCube("gs://b-cubed-eu/data_PT-30b.parquet", gproject='$GPROJECT-ID')

.. code-block:: python

    print(cube.df)

.. code-block:: text

           kingdom  kingdomkey           phylum  phylumkey              class  \
    0      Plantae           6       Charophyta    7819616       Charophyceae   
    1      Plantae           6     Tracheophyta    7707728     Polypodiopsida   
    ...
    52191  Plantae           6        Bryophyta         35          Bryopsida   
    
           classkey            order  orderkey          family  familykey  ...  \
    0           328         Charales       626       Characeae       8782  ...   
    ...
    52191       327        Pottiales       621      Pottiaceae       4671  ...   
    
          classcount  ordercount familycount genuscount distinctobservers  \
    0              1           1           1          1                 1   
    ...
    52191          4           1           1          1                 1   
    
          occurrences  mintemporaluncertainty  mincoordinateuncertaintyinmeters  \
    0               1                 2678400                            1000.0   
    ...
    52191           1                   86400                              30.0   
    
               cellCode                                           geometry  
    0      W016N32ACAAA  POLYGON ((-17 32.71875, -16.96875 32.71875, -1...  
    ...
    52191  W017N32BBCAD  POLYGON ((-17.21875 32.8125, -17.1875 32.8125,...  
    
    [52192 rows x 28 columns]

.. code-block:: python

    cube._species_richness()

.. code-block:: python

    print(cube.richness)

.. code-block:: text

                 cell  richness
    0    W016N30DDDDA         3
    1    W016N32AACAC       212
    ...
    124  W017N32BDBDD        91
    
    [125 rows x 2 columns]

.. code-block:: python

    b3cube.plot_richness(cube.richness, cube.df)

.. image:: _static/images/richness_plot.png

.. code-block:: python

    CL = griis.CheckList("$YOUR_DIRECTORY/merged_distr.txt")

.. code-block:: python

    d_s, d_c = b3cube.cumulative_species(cube, CL.species)

.. code-block:: python

    time, rate = b3cube.calculate_rate(d_c)

.. code-block:: python

    C1 = simulation.simulate_solow_costello(time, rate, vis=True)

.. code-block:: text

    Optimization terminated successfully.
         Current function value: -263.092115
         Iterations: 172
         Function evaluations: 287

.. image:: _static/images/output_9_2.png