urbanphenology
==============

Data and analysis procedures used in

Darrel Jenerette, Greg Miller, Alexander Buyantuyev, Diane E Pataki, Thomas Gillespie, Stephanie Pincetl. 
Urban vegetation dynamics and income segregation in drylands: A synthesis of seven metropolitan regions in the southwestern United States. 
Environmental Research Letters (in press)


### Data

| Variable       | Source                                                           |
|----------------|------------------------------------------------------------------|
| Growing season | MODIS 16 day 250 meter EVI (enhanced vegetation index) 2000-2010 |
| Precipitation  | Airport weather stations                                         |
| Land cover     | NLCD 2006                                                        |
| Income         | Census block income variable in shapefile format                 |


### Code Files 

File        | Description                         |
------------|-------------------------------------|
PreProc.py  | extract source data to netcdf       |
main.R      | main statistical analysis           |
util.R      | netcdf utility functions            |
phenology.R | phenology                           |
compare.R   | Compare multiple studies            |


### Create netcdf data files

### Statistical Analysis 

- phenology metrics (growing season length)
- correlation between growing season and precipitation (grouped by income and land cover)
- Compare results from several cities

### License

MIT license
 
### Author 

Greg Miller
Department of Earth Sciences 
University of California, Riverside
gmill002@gmail.com
