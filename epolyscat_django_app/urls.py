from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register("experiments", views.ExperimentViewSet, basename="experiment")
router.register("runs", views.RunViewSet, basename="run")
router.register(
    "plot-parameters", views.PlotParametersViewSet, basename="plot-parameters"
)
router.register("views", views.ViewsViewSet, basename="view")
app_name = "epolyscat_django_app"
urlpatterns = [
    path("home/", views.home, name="home"),
    path("api/plot/", views.plot, name="plot"),
    #path("api/list-inputs/", views.list_inputs, name="list_inputs"),
    #path("api/diff-inputs/", views.diff_inputs, name="diff_inputs"),
    #path("api/plotables/", views.plotables, name="plotables"),
    path("api/settings/", views.api_settings, name="settings"),
    path("api/", include(router.urls)),
]
