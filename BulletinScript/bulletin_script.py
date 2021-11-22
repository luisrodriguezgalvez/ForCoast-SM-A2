import PIL
from PIL import Image, ImageDraw,ImageFont
import os
import os.path, time
import sys, getopt
import yaml
import glob

# argv = sys.argv[1:]

# try:
#     opts, args = getopt.getopt(argv,"hy:d:s:t:",["yamlfile="])
# except getopt.GetoptError:
#     print ('bulletin_script.py -y <yamlfile>')
#     sys.exit(2)
# for opt, arg in opts:
#     if opt == '-h':
#         print ('bulletin_script.py -y <yamlfile>')
#         sys.exit()
#     elif opt in ("-y", "--yamlfile"):
#         USER_YAML_FILE = arg
#     elif opt in ("-d"):
#         datadir = arg
#     elif opt in ("-s"):
#         source = arg
#     elif opt in ("-t"):
#         target = arg

with open("../Processing/pilot7.yaml") as f:
    data = yaml.load(f, Loader=yaml.FullLoader)
    fdate=data['sdate']
    fx0=data['x0']
    fy0=data['y0']
    farea=data['experiment']

    #fdate=data['sdate']
    #fx0=data['x0']
    #fy0=data['y0']
    #farea=data['experience']

#img1 = PIL.Image.open('TS_Alarm.png')
#img2 = PIL.Image.open('FORCOAST_Logo_WhiteBack.png')

img_violin    = Image.open(figdir+'TS_violin.png')
img_risk      = Image.open(figdir+'TS_Risk.png')
img_riskchart = Image.open(figdir+'TS_Risk_chart.png')

img1 = PIL.Image.open('TS_Alarm.png')
img2 = PIL.Image.open('FORCOAST_Logo_WhiteBack.png')
img3 = PIL.Image.open('map_particles.png')
img4 = PIL.Image.open('FORCOAST_Footer_Blue.png')

riskval  = [ float(v[1]) for v in risktab]
imaxrisk = riskval.index(max(riskval))

img1Width, img1Height = img1.size
img2Width, img2Height = img2.size
img3Width, img3Height = img3.size
img4Width, img4Height = img4.size

img_violin_Width   , img_violin_Height    = img_violin.size
img_risk_Width     , img_risk_Height      = img_risk.size
img_riskchart_Width, img_riskchart_Height = img_riskchart.size
img_map_Width      , img_map_Height       = img_map.size
img_logo_Width     , img_logo_Height      = img_logo.size

img1_fixed_height = int((img1Height + img3Height) / 2)

img1_height_percent = (img1_fixed_height / float(img1.size[1]))
img1_width_size = int((float(img1.size[0]) * float(img1_height_percent)))
img1_r = img1.resize((img1_width_size, img1_fixed_height), PIL.Image.NEAREST)
# img1_r.save('img1_resized.png')


img3_fixed_height = int((img1Height + img3Height) / 2.2)

img3_height_percent = (img3_fixed_height / float(img3.size[1]))
img3_width_size = int((float(img3.size[0]) * float(img3_height_percent)))
img3_r = img3.resize((img3_width_size, img3_fixed_height), PIL.Image.NEAREST)
# img3_r.save('img3_resized.png')


img4_fixed_width = int(img1_width_size + img3_width_size)

img4_width_percent = (img4_fixed_width / float(img4.size[0]))
img4_height_size = int((float(img4.size[1]) * float(img4_width_percent)))
img4_r = img4.resize((img4_fixed_width, img4_height_size), PIL.Image.NEAREST)
# img4_r.save('footer_resized.png')


newImg = Image.new('RGBA', (img1_width_size + img3_width_size, img1_fixed_height + img2Height + img4_height_size), (255, 255, 255))

newImg.paste(img1_r, (0, img2Height))
newImg.paste(img2, (20, 20))
newImg.paste(img3_r, (img1_width_size, img2Height + int((img1_fixed_height - img3_fixed_height) / 2)))
newImg.paste(img4_r, (0, (img2Height + img1_fixed_height)))

font_path_1 = "ariali.ttf"
font_1 = ImageFont.truetype(font_path_1, 26)

font_path_2 = "arialbd.ttf"
font_2 = ImageFont.truetype(font_path_2, 40)

# print("last modified: %s" % time.ctime(os.path.getmtime(file)))
filecreated = time.ctime(os.path.getctime('TS_Alarm.png'))

draw = PIL.ImageDraw.Draw(newImg)
draw.text(((img2Height + 410, img2Height/2.6)), ('LAND POLLUTION SERVICE'), font=font_2,fill=(23,111,176,255))
draw.text((img1_width_size + img3_width_size - 600, img2Height/1), ('Bulletin generated on: ' + filecreated), font=font_1,fill=(0,0,0,255))
draw.text((img1_width_size + img3_width_size - 600, img2Height/1.4), ('Release date: ' + fdate), font=font_1,fill=(0,0,0,255))
draw.text((img1_width_size + img3_width_size - 600, img2Height/2.1), ('Discharge point: x = ' + str(fx0[0]) + ' ; ' + 'y = ' + str(fy0[0])), font=font_1,fill=(0,0,0,255))
draw.text((img1_width_size + img3_width_size - 600, img2Height/4), ('Area: ' + farea), font=font_1,fill=(0,0,0,255))

newImg.save("bulletin_with_header3.png", quality=95)




#draw = ImageDraw.Draw(newImg)
#draw.text((img1Width-400, img2Height/1)  , filecreated, font=font,fill=(0,0,0,255))
#draw.text((img1Width-400, img2Height/1.4), fdate, font=font,fill=(0,0,0,255))
#draw.text((img1Width-400, img2Height/1.9), ('x = ' + str(fx0[0]) + ' ' + 'y = ' + str(fy0[0])), font=font,fill=(0,0,0,255))
#draw.text((img1Width-400, img2Height/2.5), farea, font=font,fill=(0,0,0,255))

# newImg.save(figdir+"bulletin.png", quality=95)
'''

# draw = PIL.ImageDraw.Draw(img)

# draw.text((100, 100),"Sample Text")

# img1.paste(img2)

# img1.save("combined.png", quality=45)
'''