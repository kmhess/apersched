from astropy.coordinates import SkyCoord

# Standard calibrators:
names = ['3C138', '3C147', 'CTD93']
calibrators = [SkyCoord.from_name(name) for name in names]
