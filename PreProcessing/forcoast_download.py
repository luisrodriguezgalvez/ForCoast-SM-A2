# Download input data for A2 service module

# Input arguments required: T0 (reference date in yyyy-mm-dd) period of simulation (in days) and path where to store data
# Example: 
# - Windows: python forcoast_download.py 2021-02-01 4 c:\data
# - Linux: python forcoast_download.py 2021-02-01 4 /data

import wget
import datetime
import sys
import calendar

def download_input_variables(T0,period,datadir):

	# Need data from one day prior to T0
	delta = datetime.timedelta(days=1)
	period = datetime.timedelta(days=int(period))

	variables = ['U','V','W']

	for ii, variable in enumerate(variables):

		start_date = T0 - delta
		end_date = start_date + period + delta

		while start_date <= end_date:

			days_in_month = calendar.monthrange(int(start_date.strftime("%Y")), int(start_date.strftime("%m")))[1]

			url = 'http://46.105.99.42:8080/thredds/fileServer/PHY/' + start_date.strftime("%Y") + '/' + start_date.strftime("%m") + '/2_EFORIE_1h_' + start_date.strftime("%Y%m") + '01_' + start_date.strftime("%Y%m") + str(days_in_month) + '_grid_' + variables[ii] + '_' + start_date.strftime("%Y%m%d") + '-' + start_date.strftime("%Y%m%d") + '.nc'
			# print('http://46.105.99.42:8080/thredds/fileServer/PHY/' + start_date.strftime("%Y") + '/' + start_date.strftime("%m") + '/1_NWS_1h_' + start_date.strftime("%Y%m") + '01_' + start_date.strftime("%Y%m") + '31_grid_' + variables[ii] + '_' + start_date.strftime("%Y%m%d") + '-' + start_date.strftime("%Y%m%d") + '.nc')
			print('getting: ' + url)
			filename = wget.download(url, out=datadir)
			start_date += delta

def download_static_grids(datadir):

	files = ['2_mesh_mask.nc','mesh_mask.nc','mesh_mask_59levels.nc']

	for ii, file in enumerate(files):

		url = 'http://46.105.99.42:8080/thredds/fileServer/PHY/' + files[ii]
		print('getting: ' + url)

		filename = wget.download(url, out=datadir)

if __name__ == '__main__':

	argv = sys.argv[1:]
	T0 = argv[0]
	period = argv[1]
	datadir = argv[2]
	print('T0 = ' + str(T0))
	print('Period = ' + str(period))

	T0_datetime = datetime.datetime.strptime(str(T0),"%Y-%m-%d")
	print(T0_datetime)
	
	download_input_variables(T0_datetime,period,datadir)
	download_static_grids(datadir)
