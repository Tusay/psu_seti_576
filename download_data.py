from urllib.parse import urljoin, urlparse, parse_qs, urlencode, urlunparse
from urllib.request import urlretrieve
from mechanize import Browser
import wget
import pdb

def query_breakthrough(ra, dec):
    # parse out the url
    url = "https://breakthroughinitiatives.org/opendatasearch?project=GBT&file_type=HDF5&ra=&ra_range=&decl=&decl_range=&mjd=&mjd_range=&freq=&freq_range=&target=&perPage=&search=Search"
    parsed = urlparse(url)
    qs = parse_qs(parsed.query, keep_blank_values=True)

    # update the queries
    qs['ra'] = ra
    qs['decl'] = dec
    qs['perPage'] = 100

    # re-assemble the new url
    newqs = urlencode(qs, doseq=1)
    newurl = urlunparse([newqs if i == 4 else x for i,x in enumerate(parsed)])
    return newurl

def get_data_urls(url):
    br = Browser()
    br.set_handle_robots( False )
    br.addheaders = [('User-agent', 'Firefox')]

    # Retrieve the Google home page, saving the response
    br.open(url)
    data_urls = []
    for link in br.links(url_regex="ssl.berkeley.edu"):
        data_urls.append(link.url)
    return data_urls

def download_data_from_url(url_list, outdir):
    if type(url_list) == str:
        url_list = [url_list]
    for url in url_list:
        wget.download(url, outdir)
    return None

def download_all_data(ra, dec, outdir):
    url = query_breakthrough(ra, dec)
    data_urls = get_data_urls(url)
    download_data_from_url(data_urls, outdir)
    return None
