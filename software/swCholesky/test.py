#!/usr/bin/python

import os
import sys
import time

binFile = './sw_cholesky'
mtx_path    = '../FloridaSparse0/'
perm_path   = mtx_path + 'perm/'
perm_x86_path   = mtx_path + 'perm_x86/'
perm_log_path = 'log/log_perm/'
exec_log_path = 'log/log_exec/'

chunk = 1
costParam = 64
levelParam = -1
finalSeqNode = 4

def mkdir(path):
 
    path=path.strip()
    path=path.rstrip("/")
 
    isExists=os.path.exists(path)
 
    if not isExists:
        os.makedirs(path) 
 
        print path+' create folder'
        return True
    else:
        print path+' already exists'
        return False


def listdir(path):
    list_name = []
    for file in os.listdir(path):
        #file_path = os.path.join(path, file)
        #if os.path.isdir(file_path):
        #    listdir(file_path, list_name)
        if os.path.splitext(file)[1]=='.mtx':
            list_name.append(file)
    return list_name
def name_mtx2perm(mtx_file):
    if os.path.splitext(mtx_file)[1]=='.mtx':
        perm_file = os.path.splitext(mtx_file)[0]
        perm_file = perm_file + '.perm'
    else:
        print 'name ERROR'
    return perm_file

def gen_exec_cmd(bin_file, mtx_file, log_path, chunk, costParam, levelParam, finalSeqNode):
    perm_file = name_mtx2perm(mtx_file)
    exec_cmd = 'bsub -b -I -m 1 -p -q q_sw_expr -host_stack 1024 -sw3run ./sw3run-all -sw3runarg "-a 1" -cross_size 28000 -n 1 '
    exec_cmd += ' ' + ' -o ' + log_path + mtx_file
    exec_cmd += ' ' + '-cgsp 64'
    #exec_cmd += ' ' + './sw_cholesky'
    exec_cmd += ' ' + bin_file
    exec_cmd += ' ' + mtx_path + mtx_file
    exec_cmd += ' ' + str(chunk)
    exec_cmd += ' ' + str(costParam)
    exec_cmd += ' ' + str(levelParam)
    exec_cmd += ' ' + str(finalSeqNode)
    exec_cmd += ' ' + perm_x86_path + perm_file
    return exec_cmd


mtx_files = listdir(mtx_path)
#mtx_files = ['nd24k.mtx', 'PFlow_742.mtx']
#mtx_files = ['apache2.mtx']
#mtx_files = ['cfd2.mtx']
#mtx_files = ['bmwcra_1.mtx']
#mtx_files = ['PFlow_742.mtx']

for i in ['1']
    bin_tmp = './sw_cholesky'
    log_tmp = 'log/'
    mkdir(log_tmp)
    for mtx_file in mtx_files:
            exec_cmd = gen_exec_cmd(bin_tmp, mtx_file, log_tmp, chunk, costParam, levelParam, finalSeqNode)

            print exec_cmd
            #os.system('nohup ' + exec_cmd + ' &')
            os.system(exec_cmd)
            #time.sleep(1)
    #time.sleep(1000)










