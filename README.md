# Integrating-anthropogenic-emission-tool
This repository aims at integrating multi-source emissions as a comprehensive, standard, representative emission inventory on sectors, spatial coverage, and NMVOCs species.
Step 1. Run 'regridding.py', to create the grid proxy in 0.25 spatial degree for all data. 
Step 2. Run 'inte_all.py', to create the NETCEF files for all species with 8 sectors. 
Step 3. Run 'VOC_species_nc.py', to speciate the NMVOC emission as 20 species, which follows MOZART chemistry speciation, and finally to create the NETCEF file for every speciated emission.  
