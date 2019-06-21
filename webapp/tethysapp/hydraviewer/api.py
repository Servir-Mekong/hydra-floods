# Define your REST API endpoints here.
# In the comments below is an example.
# For more information, see:
# http://docs.tethysplatform.org/en/dev/tethys_sdk/rest_api.html
"""
from django.http import JsonResponse
from rest_framework.authentication import TokenAuthentication
from rest_framework.decorators import api_view, authentication_classes

@api_view(['GET'])
@authentication_classes((TokenAuthentication,))
def get_data(request):
    '''
    API Controller for getting data
    '''
    name = request.GET.get('name')
    data = {"name": name}
    return JsonResponse(data)
"""
