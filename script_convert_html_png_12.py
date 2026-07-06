import glob, time, os
from selenium import webdriver
from selenium.webdriver.chrome.options import Options

opts = Options()
opts.add_argument("--headless=new")
opts.add_argument("--window-size=800,800")
opts.add_argument("--hide-scrollbars")
driver = webdriver.Chrome(options=opts)

for html in sorted(glob.glob("struct_atom*_*.html")):
    png = html.replace(".html", ".png")
    url = "file:///" + os.path.abspath(html).replace("\\", "/")
    driver.get(url)
    time.sleep(3)                      # let 3Dmol/WebGL finish drawing
    driver.save_screenshot(png)
    print("wrote", png)

driver.quit()