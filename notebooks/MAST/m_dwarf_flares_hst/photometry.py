"""
Doing a little aperture photometry.

author: @arjunsavel
"""

from photutils.aperture import CircularAperture, CircularAnnulus
fluxes = []
dates = []

def flux(x, y):
    return img[x][y]# todo: check indexing!

for i, obs_id in enumerate(obs_ids):
    
    files = Observations.download_products(obs_id,
                               productType=["SCIENCE", "PREVIEW"],
                               extension="fits") 
    file = files['Local Path'][0]
    try:
        img = fits.getdata(file)
    except:
        print('!!!!!!bad!!!!!!!')
        print(i)
        continue
    
    #notice that darks, flats, skies
    header = fits.getheader(file)
    inst = header['INSTRUME']
    if inst != 'WFPC2':
        continue
    print(header['MASKCORR'], header['ATODCORR'], header['BLEVCORR'], header['BIASCORR'], header['DARKCORR'], header['FLATCORR'])
    ra_proposed = header['RA_TARG'] * u.deg # todo: grab this from the file header
    dec_proposed = header['DEC_TARG'] * u.deg # todo: grab this from the file header
    plate_scale = header['D001SCAL'] # todo: grab this from the file header from drizzle


    obs_day = Time(header['DATE-OBS'] + 'T' + header['TIME-OBS'])
    print(obs_day)

    delta_t = obs_day - J2000
    c_proposed = SkyCoord(ra=ra_proposed, dec=dec_proposed,
             obstime=J2000)

    c_target= SkyCoord(ra=df['RAJ2000'][1] * u.deg, dec=df['DEJ2000'][1] * u.deg,
                 obstime=J2000, pm_ra_cosdec=pm_ra_cosdec, pm_dec=pm_dec)

    c_target = c_target.apply_space_motion(dt=delta_t)
    sep = c_target.separation(c_proposed)
    print(sep.to(u.arcsec).value/plate_scale) # distance in arcsec!
    
    pa = position_angle(ra_proposed, df['RAJ2000'][1] * u.deg, dec_proposed, df['DEJ2000'][1] * u.deg)
    print(pa)
    
    # do photometry
    mean, median, std = sigma_clipped_stats(img, sigma=3.0)  

    fwhm=35
    daofind = DAOStarFinder(fwhm=fwhm, threshold=4.*std)  
    sources = daofind(img - median)  
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=fwhm)
    phot_table = aperture_photometry(img, apertures)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    phot_table.sort('aperture_sum')

    x = phot_table[-1]['xcenter']
    y = phot_table[-1]['ycenter']
    flux_ap = phot_table[-1]['aperture_sum']
    
    annulus_aperture = CircularAnnulus([x.value, y.value], r_in=85, r_out=95)
    aperstats = ApertureStats(img, annulus_aperture)
    bkg_mean = aperstats.mean
    

    aperture_area = apertures[0].area_overlap(img)
    total_bkg = bkg_mean * aperture_area
    
    flux_ap -= total_bkg

    fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
    im = ax.imshow(np.log10(img))
    plt.imshow(np.log10(img))
    plt.scatter(x,y)

    plt.colorbar(im)

    plt.scatter(img.shape[0]//2, img.shape[0]//2 + 24)
    plt.scatter(x,y)
    plt.show()
    
    center = (img.shape[0]//2, img.shape[0]//2 + 24) # todo: get from star-finding
    
#     aperture_size = 100

#     grid = np.where(d < aperture_size)

#     flux_ap = np.sum(flux(grid[1], grid[0]))
    fluxes += [flux_ap]
    dates += [obs_day.mjd]
    # todo: sky-subtraction
    if i == 1:
        break
    
    
    
plt.scatter(dates, fluxes)
plt.yscale('log')
plt.xlim(49750, 50000)