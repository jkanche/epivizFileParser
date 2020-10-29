from epivizfileserver import setup_app, MeasurementManager, create_fileHandler

emd_url = 'http://10.158.220.116/emd/api/v1'

mMgr = MeasurementManager()
handler = create_fileHandler()

mMgr.import_emd(emd_url, fileHandler=handler, listen=True)
print(mMgr.using_emd())

app = setup_app(mMgr)

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=8000)
