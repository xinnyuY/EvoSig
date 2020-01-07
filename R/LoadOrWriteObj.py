import gzip,pickle
import json


def load_data():
    f = gzip.open('C:\\Users\\tang\\Desktop\\Personal\\research\\coding\\data\\mnist.pkl.gz', 'rb')
    train_set, valid_set, test_set = pickle.load(f)
    return [train_set, valid_set, test_set]


def write_wat(wvh,wbv,wbh,fileName):
    f1 = open(fileName,'w')
    for x in [wvh,wbv,wbh]:
        str1 = json.dumps(x.tolist())
        f1.write("%s\n" % str1)
    f1.close()


def write_or_read_obj(path_and_name, obj=None):
    # -- Save
    if obj is not None:
        with open(path_and_name+'.pkl', 'wb') as out_f:
            pickle.dump(obj, out_f)
    # -- Read
    else:
        with open(path_and_name+'.pkl','rb') as in_f:
            obj = pickle.load(in_f)
        return obj
