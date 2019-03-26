from astropy.coordinates import SkyCoord

# Standard calibrators:
# names = ['3C138', '3C147', 'CTD93', '3C286', '3C48']
flux_names = ['3C147', '3C196']
flux_cal = [SkyCoord.from_name(name) for name in flux_names]
pol_names = ['3C138', '3C286', 'CTD93', '3C48']
pol_cal = [SkyCoord.from_name(name) for name in pol_names]

if __name__ == '__main__':
    print("Flux calibrators are:")
    for i in range(len(flux_names)):
        print('\t({}) {}: {}'.format(i, flux_names[i], flux_cal[i].to_string('hmsdms')))
    print("Polarization calibrators are:")
    for i in range(len(pol_names)):
        print('\t({}) {}: {}'.format(i, pol_names[i], pol_cal[i].to_string('hmsdms')))