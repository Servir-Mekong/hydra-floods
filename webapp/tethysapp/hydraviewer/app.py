from tethys_sdk.base import TethysAppBase, url_map_maker


class Hydraviewer(TethysAppBase):
    """
    Tethys app class for HYDRA Viewer.
    """

    name = 'HYDRAFloods Viewer'
    index = 'hydraviewer:home'
    icon = 'hydraviewer/images/hydrafloodsicon_lrgnotxt.png'
    package = 'hydraviewer'
    root_url = 'hydraviewer'
    color = '#34495e'
    description = 'Place a brief description of your app here.'
    tags = '&quot;Remote-Sensing&quot;,&quot;Floods&quot;'
    enable_feedback = False
    feedback_emails = []

    def url_maps(self):
        """
        Add controllers
        """
        UrlMap = url_map_maker(self.root_url)

        url_maps = (
            UrlMap(
                name='home',
                url='hydraviewer',
                controller='hydraviewer.controllers.home'
            ),
            UrlMap(
                name='mapviewer',
                url='hydraviewer/mapviewer',
                controller='hydraviewer.controllers.mapviewer'
            ),
            UrlMap(
                name='get_precipmap',
                url='hydraviewer/mapviewer/get_precipmap',
                controller='hydraviewer.ajax_controllers.get_precipmap'
            ),
            UrlMap(
                name='get_surfacewatermap',
                url='hydraviewer/mapviewer/get_surfacewatermap',
                controller='hydraviewer.ajax_controllers.get_surfacewatermap'
            ),
            UrlMap(
                name='download_surfacewatermap',
                url='hydraviewer/mapviewer/download_surfacewatermap',
                controller='hydraviewer.ajax_controllers.download_surfacewatermap'
            ),
            UrlMap(
                name='usecases',
                url='hydraviewer/usecases',
                controller='hydraviewer.controllers.usecases'
            ),
            UrlMap(
                name='update_historical',
                url='hydraviewer/mapviewer/update_historical',
                controller='hydraviewer.ajax_controllers.update_historical'
            ),
            UrlMap(
                name='get_latest',
                url='hydraviewer/mapviewer/get_latest',
                controller='hydraviewer.ajax_controllers.get_latest'
            ),
        )

        return url_maps
