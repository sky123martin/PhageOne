from application import app
from application.utility import setup_data

if __name__ == '__main__':
    # Threaded option to enable multiple instances for multiple user access support
    setup_data()
    app.run(threaded=True, port=5000)