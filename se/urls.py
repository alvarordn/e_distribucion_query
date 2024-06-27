from django.urls import path
from . import views

urlpatterns = [
    path('projects/', views.getProjects),
    path('project/<str:pk>', views.getProject),
    path('project/<str:pk>/solve/', views.solveProblem),
]
