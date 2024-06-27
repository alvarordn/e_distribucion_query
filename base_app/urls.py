from django.urls import path
from . import views

urlpatterns = [
    path('', views.base_view , name = 'base_view'),
    path('query/', views.endesa , name = 'endesa'),
]