import PIL
from PIL import Image, ImageDraw,ImageFont 
import os.path, time


img1 = PIL.Image.open('TS_Alarm.png')
img2 = PIL.Image.open('FORCOAST_Logo_WhiteBack.png')


img1Width, img1Height = img1.size
img2Width, img2Height = img2.size


newImg = Image.new('RGBA', (img1Width, img1Height + img2Height), (255, 255, 255))

newImg.paste(img1, (0, img2Height))
newImg.paste(img2, (15, 0))

font_path = "ariali.ttf"
font = ImageFont.truetype(font_path, 20)

# print("last modified: %s" % time.ctime(os.path.getmtime(file)))
filecreated = time.ctime(os.path.getctime('TS_Alarm.png'))

draw = PIL.ImageDraw.Draw(newImg)
draw.text((img1Width-400, img2Height/2),filecreated, font=font,fill=(0,0,0,255))

newImg.save("bulletin_with_header2.png", quality=95)






# draw = PIL.ImageDraw.Draw(img)

# draw.text((100, 100),"Sample Text")

# img1.paste(img2)

# img1.save("combined.png", quality=45)

        



        
        
        
