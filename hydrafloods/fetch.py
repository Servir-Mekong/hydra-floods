import netrc
import ftplib
import datetime

def viirs(date,h,v,outdir='./',creds=None):
    """Function to download VIIRS NRT data for specified time and tile

    Args:
        date (datetime.datetime): Datetime object specifying which date the data of interest was acquired.
        h (int): horizontal tile grid to fetch
        v (int): vertical tile grid to fetch
        outdir (str, optional): out directory to dump retrieved data to
        default = './' or current working directory
        creds (str, optional): path to .netrc file with NASA EarthData login in credentials
        default = None

    Returns:
        None
    """

    if outdir[-1] != '/':
        outdir = outdir+'/'

    acct = netrc.netrc(creds)
    usr,_,pswrd = acct.hosts['https://urs.earthdata.nasa.gov']

    basename = 'VNP09GA_NRT.A{0}{1:03d}.h{2:02d}v{3:02d}.001.h5'

    yr = date.year
    dt = (date-datetime.datetime(yr,1,1)).days + 1

    url = 'nrt3.modaps.eosdis.nasa.gov'
    directory = '/allData/5000/VNP09GA_NRT/Recent/'

    ftp = ftplib.FTP(url)
    ftp.login(usr,pswrd)

    ftp.cwd(directory)

    filename = basename.format(yr,dt,h,v)
    with open(outdir + filename, 'wb') as f:
       ftp.retrbinary('RETR ' + filename, f.write)

    ftp.close()

    return

def atms(date,outdir='./'creds=None):
    if outdir[-1] != '/':
        outdir = outdir+'/'
        
    
    
    return

if __name__ == "__main__":
    outdir = '/home/ubuntu/hydra/data/'
    creds = '/home/ubuntu/hydra/earthdata.netrc'
    h = 27
    v = 7
    now = datetime.datetime.now() - datetime.timedelta(1)
    fetchNrtViirs(now,h,v,outdir,creds)
