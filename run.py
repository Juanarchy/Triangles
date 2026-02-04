from main import execute

options={"nxny":[100,10],
         "folder":"test_rose",
         "Tmax":180,
         "tol_dry":0.000001,
         "g":9.81,
         "manning":0,
         "CFL":0.25,
         "dt_save":1,
         "bconds":["soft","soft","periodic","periodic"], #left,right,top,bottom (west, east, north, south)
         "divisions":0,
         "triangles":"equilateral",
         "nestings":0,
         "test":"dambreak_rose"}

execute(options)