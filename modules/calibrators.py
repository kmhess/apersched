from astropy.coordinates import SkyCoord

# Standard calibrators:
names = ['3C138', '3C147', 'CTD93', '3C286']
calibrators = [SkyCoord.from_name(name) for name in names]

if __name__ == '__main__':
    for i in range(len(names)):
        print('({}) {}: {}'.format(i,names[i],calibrators[i]))