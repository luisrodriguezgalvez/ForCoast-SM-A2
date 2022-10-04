#(C) Copyright FORCOAST H2020 project under Grant No. 870465. All rights reserved

import PIL
from PIL import Image, ImageDraw,ImageFont 
import os
import os.path, time
import sys, getopt
import yaml
import csv
import glob
import cv2

def resize_width(image, width, height, new_width):
    new_height = int((new_width/width)*height)
    image_resize = image.resize((new_width, new_height), PIL.Image.NEAREST)
    return image_resize, new_width, new_height

def make_video(frame_folder, bull_width, bull_height):
    images = sorted(glob.glob(f"{frame_folder}/bulletin_*.png"))
    for i, image in enumerate(images):
            ResizeBulletin = Image.new('RGBA', (1920, 1080), (0, 0, 0))
            Bulletin = Image.open(image)
            Bulletin_resize, Bulletin_resize_Width, Bulletin_resize_Height = resize_width(Bulletin, bull_width, bull_height, 1920)
            ResizeBulletin.paste(Bulletin_resize, (0, 540 - (int(Bulletin_resize_Height/2))))
            if i < 10:
                ResizeBulletin.save("{}/0{}_resize.png".format(frame_folder, i), quality = 95)
            else:
                ResizeBulletin.save("{}/{}_resize.png".format(frame_folder, i), quality = 95)
            
    imagesResize = [cv2.imread(imageResize) for imageResize in sorted(glob.glob(f"{frame_folder}/*_resize.png"))]
    if len(images) == 0:
        print("No frames for video")
    else:
        fourccTelegram = cv2.VideoWriter_fourcc(*'mp4v')
        TelegramVideo = cv2.VideoWriter(frame_folder+"bulletin.mp4", fourccTelegram, 0.7, (1920, 1080))
        print("Telegram Video Writer Initiatied")
        for imageResize in imagesResize:
            TelegramVideo.write(imageResize)
        print("Telegram Video generated")

        fourccWeb = cv2.VideoWriter_fourcc(*'vp80')
        Webvideo = cv2.VideoWriter(frame_folder+"bulletin.webm", fourccWeb, 0.7, (1920, 1080))
        print("Web Video Writer Initiatied")
        for imageResize in imagesResize:
            Webvideo.write(imageResize)
        print("Web Video generated")

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


# One Bulletin

img_violin    = Image.open(figdir+'TS_violin.png')
img_risk      = Image.open(figdir+'TS_Risk.png')
img_riskchart = Image.open(figdir+'TS_Risk_chart.png')

risktab=[]
# Here we need to select the map for the worst case
with open(figdir+'Risk.csv', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader: # each row is a list
        risktab.append(row)

#print(risktab)

riskval  = [ float(v[1]) for v in risktab]
imaxrisk = riskval.index(max(riskval))

img_map = Image.open(figdir+'AllTracks_Alarm_'+'%03d'%(imaxrisk)+'.png')
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

# One Bulletin for each timeframe

for ti,t in enumerate(risktab):
    img_violin    = Image.open(figdir+'TS_violin_'+'%03d'%(ti)+'.png')
    img_risk      = Image.open(figdir+'TS_Risk_'+'%03d'%(ti)+'.png')
    img_riskchart = Image.open(figdir+'TS_Risk_chart.png')

    img_map = Image.open(figdir+'AllTracks_Alarm_'+'%03d'%(ti)+'.png')
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

    newImg.save(figdir+'bulletin_'+'%03d'%(ti)+'.png', quality = 95)

make_video(figdir, (2 * margin + img_violin_Width + img_map_Width), (margin + img_logo_new_Height + img_violin_Height + img_risk_Height + img_footer_Height))



