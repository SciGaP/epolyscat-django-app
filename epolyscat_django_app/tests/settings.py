
import os

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SECRET_KEY = "abc123"
INSTALLED_APPS = [
    'django.contrib.auth',
    'django.contrib.contenttypes',
    "epolyscat_django_app.apps.TrecxDjangoAppConfig",
]
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASEDIR, 'db.sqlite3'),
    }
}
