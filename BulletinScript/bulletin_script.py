import PIL
from PIL import Image, ImageDraw,ImageFont 
import os
import os.path, time
import sys, getopt
import yaml
import csv
import glob

rdate = 0

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hy:s:c:t:T:k:d:",["yamlfile="])
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
    elif opt in ("-T"):
        rdate = arg
    elif opt in ("-k"):
        targetcount = int(arg)
    elif opt in ("-d"): # AC 10112021 shouldn't be needed ? 
        datadir = arg
USER_YAML_FILE = '../usr/' + USER_YAML_FILE + '/config/config.yaml'
print ('yaml file is', USER_YAML_FILE)

with open(USER_YAML_FILE) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)
    figdir        = '../usr/'+data['username'] +'/output/target_' + str(targetcount)+'_source_' + str(sourcecount)+'/'
    
    if rdate == 0:
        rdate = data['sdate']
        
    farea=data['experiment']
    #fx0=data['x0']
    #fy0=data['y0']


img_violin    = Image.open(figdir+'TS_violin.png')
img_risk      = Image.open(figdir+'TS_Risk.png')
img_riskchart = Image.open(figdir+'TS_Risk_chart.png')

risktab=[]
# Here we need to select the maps for the worst case
with open(figdir+'Risk.csv', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader: # each row is a list
        risktab.append(row)

print(risktab)

riskval  = [ float(v[1]) for v in risktab]
imaxrisk = riskval.index(max(riskval))

img_map = Image.open(figdir+'AllTracks_Alarm'+str(imaxrisk)+'.png')
img_logo = Image.open('./FORCOAST_Logo_WhiteBack.png')
img_footer = Image.open('./FORCOAST_Footer_Blue.png')

img_violin_Width   , img_violin_Height    = img_violin.size
img_risk_Width     , img_risk_Height      = img_risk.size
img_riskchart_Width, img_riskchart_Height = img_riskchart.size
img_map_Width      , img_map_Height       = img_map.size
img_logo_Width     , img_logo_Height      = img_logo.size
img_footer_Width   , img_footer_Height    = img_footer.size

margin = 25

# Resize the logo

img_logo_new_Height = 220

img_logo_height_percent = (img_logo_new_Height / float(img_logo.size[1]))
img_logo_new_Width = int((float(img_logo.size[0]) * float(img_logo_height_percent)))
img_logo_new = img_logo.resize((img_logo_new_Width, img_logo_new_Height), PIL.Image.NEAREST)

# Generate the new combined image

newImg = Image.new('RGBA', (2 * margin + img_violin_Width + img_map_Width, margin + img_logo_new_Height + img_violin_Height + img_risk_Height + img_footer_Height), (255, 255, 255))

newImg.paste(img_logo_new , (margin                                     , margin))
newImg.paste(img_violin   , (margin                                      , margin + img_logo_new_Height))
newImg.paste(img_map      , (margin + img_violin_Width                       , margin + img_logo_new_Height))
newImg.paste(img_risk     , (margin                                      , margin + img_logo_new_Height + img_violin_Height))
newImg.paste(img_riskchart, (margin + img_violin_Width                       , margin + img_logo_new_Height + img_violin_Height))
newImg.paste(img_footer   , (margin + int((img_violin_Width + img_map_Width) / 2 - (img_footer_Width /2)), margin + img_logo_new_Height + img_violin_Height + img_risk_Height))

font_path_1 = "ariali.ttf"
font_1 = ImageFont.truetype(font_path_1, 36)

font_path_2 = "ariali.ttf"
font_2 = ImageFont.truetype(font_path_2, 26)

font_path_3 = "arialbd.ttf"
font_3 = ImageFont.truetype(font_path_3, 60)

# print("last modified: %s" % time.ctime(os.path.getmtime(file)))
filecreated = time.ctime(os.path.getctime(figdir+'TS_violin.png'))

draw = PIL.ImageDraw.Draw(newImg)
draw.text(((img_logo_new_Width + 410, img_logo_new_Height / 2.6)), ('LAND POLLUTION SERVICE'), font=font_3,fill=(23,111,176,255))
draw.text((img_violin_Width + img_map_Width - 750, img_logo_new_Height / 1.1), ('Bulletin generated on: ' + filecreated), font=font_2,fill=(0,0,0,255))
draw.text((img_violin_Width + img_map_Width - 750, img_logo_new_Height / 2.1), ('Release date: ' + rdate), font=font_1,fill=(0,0,0,255))
draw.text((img_violin_Width + img_map_Width - 750, img_logo_new_Height / 4), ('Area: ' + farea), font=font_1,fill=(0,0,0,255))
# draw.text((img_violin_Width + img_map_Width - 600, img_logo_new_Height / 5), ('x = ' + str(fx0[0]) + ' ' + 'y = ' + str(fy0[0])), font=font_1,fill=(0,0,0,255))


newImg.save(figdir+"bulletin.png", quality = 95)

'''
# draw = PIL.ImageDraw.Draw(img)
# draw.text((100, 100),"Sample Text")
# img1.paste(img2)
# img1.save("combined.png", quality=45)
'''
