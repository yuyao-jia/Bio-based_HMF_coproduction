def convert_listdata_to_csv(fname='data',listdata=None):
    with open(f"{fname}.csv", 'w+') as tf:
        for idx, _arr in enumerate(listdata):
            tf.write(f"arr_{idx},{','.join(str(val) for val in _arr)}\n")