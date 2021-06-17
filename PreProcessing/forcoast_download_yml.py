# Download input data for A2 service module

# Input arguments required: T0 (reference date in yyyy-mm-dd) period of simulation (in days) and path where to store data
# Example: 
# - Windows: python forcoast_download.py 2021-02-01 4 c:\data
# - Linux: python forcoast_download.py 2021-02-01 4 /data

import wget
import datetime
import sys
import calendar
import yaml
import ftplib
from urllib.parse import urlparse

def download_files(SM,pilot,T0,period,download_dict):

	# Loop over all pilot areas in yaml file
	for ii in range(len(download_dict)):

		# Only download data for relevant Service Module and pilot
		if download_dict[ii][ii]['pilot'] == pilot and download_dict[ii][ii]['service_module'] == SM:

			# Check if output data should be downloaded
			if download_dict[ii][ii]['download_data']:

				# Set time information
				period = datetime.timedelta(days=int(period))
				delta = datetime.timedelta(days=1)
				# Some data required 1 day before and after actual time window to be downloaded
				start_date = T0 - datetime.timedelta(days=download_dict[ii][ii]['time_offset'])
				end_date = start_date + period + datetime.timedelta(days=download_dict[ii][ii]['time_offset'])

				duration = end_date - start_date

				for tt in range(0,duration.days+1):
					print(tt)
					for jj in range(len(download_dict[ii][ii]['datafiles'])):
					
						url = download_dict[ii][ii]['datafiles'][jj]

						# day of year
						day_of_year = (start_date - datetime.datetime(start_date.year, 1, 1)).days + 1 

						# Replace all date and time info in URL's
						url = url.replace('(YYYY)',start_date.strftime("%Y"))
						url = url.replace('(YYYYmm)',start_date.strftime("%Y%m"))
						url = url.replace('(YYYYmmdd)',start_date.strftime("%Y%m%d"))
						url = url.replace('(YYYYmmddHH)',start_date.strftime("%Y%m%d%H"))
						url = url.replace('(mm)',start_date.strftime("%m"))
						url = url.replace('(ddd)',str(day_of_year))
						url = url.replace('(month_length)',str(calendar.monthrange(start_date.year, start_date.month)[1]))
						
						if download_dict[ii][ii]['method'] == "http":

							print('getting: ' + url)
							filename = wget.download(url, out=download_dict[ii][ii]['outpath'])

						# Switching to ftplib since wget doesnt seem to support username and password properly
						elif download_dict[ii][ii]['method'] == "ftp":

							url_split = url.split("/")[-15:]

							ftp_ip = url_split[2]
							ftp_file = url_split[len(url_split)-1]

							ftp_folder = ''
							for xx in range(3, len(url_split)-1):
								ftp_folder = ftp_folder + '/' + url_split[xx]

							ftp = ftplib.FTP(ftp_ip) 
							ftp.login(download_dict[ii][ii]['username'], download_dict[ii][ii]['password']) 
							ftp.cwd(ftp_folder)
							ftp.retrbinary("RETR " + ftp_file, open((download_dict[ii][ii]['outpath'] + '/' + ftp_file), 'wb').write)
							ftp.quit()

					start_date += delta

			# Check if grid definition should be downloaded as well
			if download_dict[ii][ii]['download_grid']:
				for jj in range(len(download_dict[ii][ii]['gridfiles'])):
					url = download_dict[ii][ii]['gridfiles'][jj]
					
					print('getting: ' + url)
					filename = wget.download(url, out=download_dict[ii][ii]['outpath'])

if __name__ == '__main__':

	argv = sys.argv[1:]

	# Starttime and period provided as input arguments
	# TODO: set to now if no input is provided
	SM = argv[0]
	pilot = argv[1]
	T0 = argv[2]
	period = argv[3]

	print('T0 = ' + str(T0))
	print('Period = ' + str(period))

    # Get further info on data to download from yaml file
	download_yaml = open("forcoast_download.yaml")
	download_dict = yaml.load(download_yaml, Loader=yaml.FullLoader)

    # Convert T0 to datetime
	T0_datetime = datetime.datetime.strptime(str(T0),"%Y-%m-%d")
	
	download_files(SM,pilot,T0_datetime,period,download_dict)
