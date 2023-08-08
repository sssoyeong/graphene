from django.conf.urls import url
# from ml import views
# from predict_prdf_NN.ml import views
from . import views

urlpatterns = [
    url(r'^$', views.IndexView.as_view(), name='index'),
]
