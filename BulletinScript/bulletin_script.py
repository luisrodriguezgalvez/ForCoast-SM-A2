#import PIL
from PIL import Image, ImageDraw,ImageFont 
import os.path, time
import sys, getopt
import yaml
import csv

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hy:s:c:t:k:d:",["yamlfile="])
except getopt.GetoptError:
    print(opts)
    print ('USAGE : SM-A2-Postprocess.py -y <username>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('USAGE : SM-A2-Postprocess.py -y <username>')
        sys.exit()
    elif opt in ("-y", "--yamlfile"):
        USER_YAML_FILE = arg
    elif opt in ("-s"):
        # AC 10112021 - Not sure this is a safe parsing method .. Same below.
        exec('source = '+arg)
    elif opt in ("-c"):
        sourcecount = int(arg)
    elif opt in ("-t"):
        exec('target = '+arg)
    elif opt in ("-k"):
        targetcount = int(arg)
    elif opt in ("-d"): # AC 10112021 shouldn't be needed ? 
        datadir = arg
USER_YAML_FILE = '../usr/' + USER_YAML_FILE + '/config/config.yaml'
print ('yaml file is', USER_YAML_FILE)

with open(USER_YAML_FILE) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)
    figdir        = '../usr/'+data['username'] +'/output/target_' + str(targetcount)+'_source_' + str(sourcecount)+'/'

    #fdate=data['sdate']
    #fx0=data['x0']
    #fy0=data['y0']
    #farea=data['experience']

#img1 = PIL.Image.open('TS_Alarm.png')
#img2 = PIL.Image.open('FORCOAST_Logo_WhiteBack.png')

img_violin    = Image.open(figdir+'TS_violin.png')
img_risk      = Image.open(figdir+'TS_Risk.png')
img_riskchart = Image.open(figdir+'TS_Risk_chart.png')

risktab=[]
# Here we need to select the maps for the worst case
with open(figdir+'Risk.csv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader: # each row is a list
        risktab.append(row)

riskval  = [ float(v[1]) for v in risktab]
imaxrisk = riskval.index(max(riskval))

img_map = Image.open(figdir+'AllTracks_Alarm'+str(imaxrisk)+'.png')
img_logo = Image.open('./FORCOAST_Logo_WhiteBack.png')

img_violin_Width   , img_violin_Height    = img_violin.size
img_risk_Width     , img_risk_Height      = img_risk.size
img_riskchart_Width, img_riskchart_Height = img_riskchart.size
img_map_Width      , img_map_Height       = img_map.size
img_logo_Width     , img_logo_Height      = img_logo.size

newImg = Image.new('RGBA', (img_violin_Width + img_map_Width, img_violin_Height + img_risk_Height), (255, 255, 255))

newImg.paste(img_violin   , (0                , 0))
newImg.paste(img_map      , (img_violin_Width , 0))
newImg.paste(img_risk     , (0                , img_violin_Height))
newImg.paste(img_riskchart, (img_violin_Width , img_violin_Height))

font_path = "ariali.ttf"
font = ImageFont.truetype(font_path, 20)

# print("last modified: %s" % time.ctime(os.path.getmtime(file)))
filecreated = time.ctime(os.path.getctime(figdir+'TS_violin.png'))

#draw = ImageDraw.Draw(newImg)
#draw.text((img1Width-400, img2Height/1)  , filecreated, font=font,fill=(0,0,0,255))
#draw.text((img1Width-400, img2Height/1.4), fdate, font=font,fill=(0,0,0,255))
#draw.text((img1Width-400, img2Height/1.9), ('x = ' + str(fx0[0]) + ' ' + 'y = ' + str(fy0[0])), font=font,fill=(0,0,0,255))
#draw.text((img1Width-400, img2Height/2.5), farea, font=font,fill=(0,0,0,255))

newImg.save(figdir+"bulletin.png", quality=95)
'''

# draw = PIL.ImageDraw.Draw(img)

# draw.text((100, 100),"Sample Text")

# img1.paste(img2)

# img1.save("combined.png", quality=45)
'''