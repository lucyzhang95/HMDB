# import re
#
# import requests
# from bs4 import BeautifulSoup


def get_release(self):
    # # 07/08/2024: version is 5.0 and since 2021-11-17
    # url = "https://hmdb.ca/downloads"
    # rep = requests.get(url)
    # content = rep.content
    #
    # parser = BeautifulSoup(content, "html.parser")

    # version_tag = parser.find("a", {"href": "#version_5"}).text
    # version = re.search(r"\((.+)\)", version_tag).group(1) if version_tag else "Version not found"
    # return str(version)
    return "5.0"
