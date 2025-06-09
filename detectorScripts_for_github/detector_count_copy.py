#Code used with NERSC Perlmutter to extract the number of detectors used to observe a GW event
#900.fits is a locally-accessed test FITS skymap

from astropy.io import fits
address = r"/Users/mck/Desktop/palmese_research/900.fits"

def detectorsCount(address):
    with fits.open(address) as hduList:
        metadata = hduList[1].header
        detectors = metadata["INSTRUME"]
        return detectors, len(detectors.split(","))
        
print(detectorsCount(address))