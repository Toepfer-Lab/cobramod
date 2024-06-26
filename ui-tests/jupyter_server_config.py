from jupyterlab.galata import configure_jupyter_server

configure_jupyter_server(c)

from tempfile import mkdtemp

c.ServerApp.port = 8888
c.ServerApp.token = ""
c.ServerApp.password = ""
c.ServerApp.disable_check_xsrf = True
c.ServerApp.root_dir = mkdtemp(prefix="galata-test-")
c.ServerApp.open_browser = False
c.LabApp.expose_app_in_browser = True
