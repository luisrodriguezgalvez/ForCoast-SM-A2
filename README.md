# ForCoast-SM-A2

### Description

This service provides a measure of the probability for farming sites to be affected by harmfull land discharge. A source and farming location  must be selected and a simulation start date must be given. The pollutants dispersion is forecasted for up to 10 days. The result is a bulletin showing the pollutants track, details on pollutants in the selected farming location and a risk indication for the farm (see example bulletin).

### How to run

* Containerize contents in docker
* Run the command Docker run forcoast/forcoast-sm-a1 &lt;pilot> &lt;date> &lt;period> &lt;mode> &lt;config> &lt;sources> &lt;targets> &lt;datadir> &lt;Telegram token> &lt;Telegram chat_id>
  * Where <pilot> is either "galway", "venice" or "eforie"
  * Where period is in days, can go as far as your data range. Usually 3 is given
  * Mode is how the sources and targets input are given
  * Config is the location of the config file
  * Sources is the pollutant source
  * Targets is a rectangle for your farm location
  * Telegram bot is used for sendingh the bulletins through messaging services
* Example of use: Docker run forcoast/forcoast-sm-a2 eforie 2022-07-09 3 github https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/Eforie_case_1/config.yaml https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/Eforie_case_1/sources.txt https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/Eforie_case_1/targets.txt /usr/src/app/data/ 5267228188:AAGx60FtWgHkScBb3ISFL1dp6Oq_9z9z0rw -1001780197306

### Licence

Licensed under the Apache License, Version 2.0
