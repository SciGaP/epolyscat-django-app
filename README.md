

# ePolyscat Django App

## Getting Started

Follow the instructions for installing the [Airavata Django Portal](https://github.com/apache/airavata-django-portal).

With the Airavata Django Portal virtual environment activated, clone this repo and install it into the portal's virtual environment:

```bash
git clone --recursive https://github.com/SciGaP/epolyscat-django-app.git
cd epolyscat-django-app
export NODE_OPTIONS=--openssl-legacy-provider
pip install -e .
````

Start (or restart) the Django Portal server.

---

## Frontend Development

To build or develop the Vue.js frontend code, follow these steps:

1. Install **Node 16**
2. Install **Yarn**
3. In `./epolyscat_django_app/`, run:

```bash
cd ./epolyscat_django_app
yarn
yarn build
```

Or, to start the dev server for local development:

```bash
cd ./epolyscat_django_app
yarn
yarn serve
```

### Running the frontend with the Django Portal

The Vue development server and the Django Portal are separate processes.
Running `yarn serve` starts the Vue development server for frontend development
and hot reloading; it does not start the Django Portal or the Airavata backend.

To test ePolyScat as an integrated application:

1. Complete the Getting Started instructions and install ePolyScat into the
   Django Portal's Python virtual environment.
2. Configure the Portal with valid Airavata API and authentication settings.
3. Start or restart the Django Portal and ensure its required services are
   available.
4. In a separate terminal, start the Vue development server using the
   Frontend Development instructions above.
5. Open ePolyScat through the application entry in the Django Portal. Do not
   assume that opening the Vue development server's root URL provides the
   complete integrated application.

#### Troubleshooting

- **`AiravataAPI is not defined`:** Check that ePolyScat is being loaded
  through the Portal environment and that the Airavata API JavaScript library
  is available before the application code executes.
- **API requests return 404 or pages remain loading:** Check that the Portal
  and the required backend services are running and that the requested API
  routes are available.
- **The frontend build succeeds but the application does not work:** A
  successful Vue build verifies frontend compilation; it does not verify
  authentication, Portal integration, or backend requests.

---

## Creating DB Migrations

```bash
django-admin makemigrations --pythonpath . --settings epolyscat_django_app.tests.settings epolyscat_django_app
```
