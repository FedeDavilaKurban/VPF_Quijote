file_parameters_vpf = dict(
    {   
        "ns":int(10000),
        "rbin":int(13),
        "rmin":3.,
        "rmax":int(25),
        "njk":int(10),
        "space":'zspace',
        "filedir":'../data/output/'
    })

file_parameters_2pcf = dict(
    {   
        "nbins_s":9,  #9, 12, 29
        "nbins_m":10,  #10, 30
        "rmin":10,     #3, 10
        "rmax":150,    #25, 150
        "space":'zspace',
        "filedir":'../data/output/2pcf/',
        "axis":0
    })