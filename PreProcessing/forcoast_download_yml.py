# Download input data for service modules and pilot areas

# Input arguments:
# 		- required: T0 (reference date in yyyy-mm-dd)  and path where to store data
#		- optional: period of simulation in days (defaults to zero)
#       - optional: service module ID (defaults to all)
#       - optional: pilot area ID (defaults to all)
# Examples: 
# - python forcoast_download_yml.py -T 2021-06-24 -p 2 -s A2 -a eforie --> download 2 days of data for A2 service module for Eforie pilot area
# - python forcoast_download_yml.py -T 2021-06-24 -p 2 --> download 2 days of data for all service modules and pilot areas enabled in yml
# - python forcoast_download_yml.py -T 2021-06-24 -p 2 -s A2 --> download 2 days of data for A2 service modules and pilot areas enabled in yml

import wget
import datetime
import sys
import calendar
import yaml
import ftplib
from urllib.parse import urlparse
import argparse
import os
import time
from pathlib import Path
# from urllib import request
import urllib.request
import ssl

def download_files(SM,pilot,T0,period,download_dict,datadir):

	# Loop over all pilot areas in yaml file
	for ii in range(len(download_dict)):

		# Only download data for relevant Service Module and pilot
		# Cases: 
		# - pilot area of service module are specified in input
		# - only pilot area specified, so download input for all service modules
		# - only service module specified, so download input all for pilot areas
		# - neither service module nor pilot area specified, so download all data
		# if pilot == "all" and SM == "all":
		if (download_dict[ii][ii]['pilot'] == pilot and download_dict[ii][ii]['service_module'] == SM) or \
		   (download_dict[ii][ii]['pilot'] == pilot and SM == "all") or \
		   (pilot == "all" and download_dict[ii][ii]['service_module'] == SM) or \
		   (pilot == "all" and SM == "all"):

		   	# Check if output data should be downloaded
			if download_dict[ii][ii]['download_data']:

				# Set time information
				delta = datetime.timedelta(days=1)
				# Some data required 1 day before and after actual time window to be downloaded
				start_date = T0 - datetime.timedelta(days=download_dict[ii][ii]['time_offset'])
				end_date = start_date + datetime.timedelta(days=int(period)) + datetime.timedelta(days=download_dict[ii][ii]['time_offset'])

				duration = end_date - start_date

				for tt in range(0,duration.days+1):
					for jj in range(len(download_dict[ii][ii]['datafiles'])):
					
						url = download_dict[ii][ii]['datafiles'][jj]

						# day of year
						day_of_year = (start_date - datetime.datetime(start_date.year, 1, 1)).days + 1 

						# Replace all date and time info in URL's
						url = url.replace('(YYYY)',start_date.strftime("%Y"))
						url = url.replace('(YYYYmm)',start_date.strftime("%Y%m"))
						url = url.replace('(YYYYmmdd)',start_date.strftime("%Y%m%d"))
						url = url.replace('(YYYY-mm-dd)',start_date.strftime("%Y-%m-%d"))
						url = url.replace('(YYYYmmdd_T0)',T0.strftime("%Y%m%d"))
						url = url.replace('(YYYYmmddHH)',start_date.strftime("%Y%m%d%H"))
						url = url.replace('(mm)',start_date.strftime("%m"))
						url = url.replace('(ddd)',str(day_of_year))
						url = url.replace('(month_length)',str(calendar.monthrange(start_date.year, start_date.month)[1]))

						# Replace data directory in yaml file with one from command line, if provided
						if datadir != "Using outpath from yml file":
							download_dict[ii][ii]['outpath'] =  datadir

						if download_dict[ii][ii]['method'] == "urllib":

							print('urllib: getting: ' + url)
							
							# Use to avoid problems with ssl verification (https://stackoverflow.com/questions/33770129/how-do-i-disable-the-ssl-check-in-python-3-x)
							ctx = ssl.create_default_context()
							ctx.check_hostname = False
							ctx.verify_mode = ssl.CERT_NONE

							filename = url.split("/")[-1:]
							output_filename = Path(download_dict[ii][ii]['outpath'] + '\\' + filename[0].replace(':',''))

							with urllib.request.urlopen(url, context=ctx) as u, \
								open(output_filename, 'wb') as f:
								f.write(u.read())
						
						if download_dict[ii][ii]['method'] == "http":

							print('\n wget: getting: ' + url)

							try:
								filename = wget.download(url, out=download_dict[ii][ii]['outpath'])
							except:
								print('File is not available, check available dates at source')

						# Switching to ftplib since wget doesnt seem to support username and password properly
						elif download_dict[ii][ii]['method'] == "ftp":

							try:

								print('\n ftp: getting: ' + url)

								url_split = url.split("/")[-15:]

								ftp_ip = url_split[2]
								ftp_file = url_split[len(url_split)-1]

								ftp_folder = ''
								for xx in range(3, len(url_split)-1):
									ftp_folder = ftp_folder + '/' + url_split[xx]

								ftp = ftplib.FTP(ftp_ip) 
								ftp.login(download_dict[ii][ii]['username'], download_dict[ii][ii]['password']) 
								ftp.cwd(ftp_folder)

								outputfile = download_dict[ii][ii]['outpath'] + '/' + ftp_file

								# Only download if file exists
								if ftp_file in ftp.nlst():
									ftp.retrbinary("RETR " + ftp_file, open(outputfile, 'wb').write)
								else:
									print('File is not available, check available dates at source')

								ftp.quit()

							except:
								
								print("Problem downloading data for Service Module " + SM + "and pilot area" + pilot)

					start_date += delta

			# Check if grid definition should be downloaded as well
			if download_dict[ii][ii]['download_grid']:
				for jj in range(len(download_dict[ii][ii]['gridfiles'])):
					url = download_dict[ii][ii]['gridfiles'][jj]

					# Replace data directory in yaml file with one from command line, if provided
					if datadir != "Using outpath from yml file":
						download_dict[ii][ii]['outpath'] =  datadir
					
					print('getting: ' + url)
					filename = wget.download(url, out=download_dict[ii][ii]['outpath'])

if __name__ == '__main__':

	# Get input from command line arguments

	parser = argparse.ArgumentParser(description = "Description for my parser")
	parser.add_argument("-T", "--time", help = "Set reference time", required = True, default = "")
	parser.add_argument("-s", "--service", help = "Service Module ID", required = False, default = "all")
	parser.add_argument("-a", "--pilot", help = "Pilot area ID", required = False, default = "all")
	parser.add_argument("-d", "--datadir", help = "Data directory in which to download data", required = False, default = "Using outpath from yml file")
	parser.add_argument("-p", "--period", help = "Set period for which to download data", required = False, default = "0")

	argument = parser.parse_args()
	status = False

	# Get input from yaml file
	download_yaml = open("forcoast_download.yaml")
	download_dict = yaml.load(download_yaml, Loader=yaml.FullLoader)

	if argument.time:
		T0 = argument.time
		print('T0 = ' + T0)
		# Convert T0 to datetime
		T0_datetime = datetime.datetime.strptime(str(T0),"%Y-%m-%d")
	if argument.period:
		period = argument.period
		print('Period = ' + period)
	if argument.service:
		SM = argument.service
		print('Service Module = ' + SM)
	if argument.pilot:
		pilot = argument.pilot
		print('Pilot area = ' + pilot)
	if argument.datadir:
		datadir = argument.datadir
		print('Data directory = ' + datadir)
		# Check if directory exists, and make it otherwise:
		isExist = os.path.exists(datadir)
		if not isExist and datadir != "Using outpath from yml file":
			os.makedirs(datadir)
			print("Data directory created: " + datadir)

	download_files(SM,pilot,T0_datetime,period,download_dict,datadir)






  
  # Create a new directory because it does not exist 

