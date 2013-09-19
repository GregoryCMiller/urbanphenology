:: Windows batch script to process all study areas
:: turn off sleep and standby modes
:: run each area, pipe output to <x>log.txt 
powercfg -change -hibernate-timeout-ac 0
powercfg -change -standby-timeout-ac 0

python create_netcdf.py config\\common.yaml config\\lbb.yaml > lbblog.txt 2>&1
python create_netcdf.py config\\common.yaml config\\phx.yaml > phxlog.txt 2>&1
python create_netcdf.py config\\common.yaml config\\aus.yaml > auslog.txt 2>&1
python create_netcdf.py config\\common.yaml config\\las.yaml > laslog.txt 2>&1
python create_netcdf.py config\\common.yaml config\\lax.yaml > laxlog.txt 2>&1
python create_netcdf.py config\\common.yaml config\\ie.yaml  > ielog.txt  2>&1
python create_netcdf.py config\\common.yaml config\\abq.yaml > abqlog.txt 2>&1
