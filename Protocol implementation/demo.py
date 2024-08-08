from builtins import str
import hashlib
import hmac
import random
import timeit
import time
import socket
import struct
import random
import datetime
from Components import Complex
from Components import EC, HalfComplex
from Components import Math
from Components.BonehFranklin import BonehFranklin_stream, BonehFranklin_block, point_to_byte, byte_to_point
from Components.ECPoint import ECPoint
from Components.Fp2Point import Fp2Point
from Components.Fp2Element import Fp2Element
from Components.TatePairing import TatePairing
from Components.HalfComplex import HalfComplex

from ecdsa import SigningKey, VerifyingKey, NIST384p



LOG = True
VERBOSE = True
logs_time = open("time.log","w")
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


def bilinear_map(gammaQv, Qr, ec):
    return TatePairing.computeF(TatePairing, gammaQv, Qr, ec)

def point_from_string(s):
    s = s.split(" ")
    return Fp2Point(int(s[0]), int(s[1]))

def registry_rnd_car(lp,gamma,car_id = None):
    start = timeit.default_timer()
    Qv = ec.at(lp)
    gQv = ec.mulJ(Qv, gamma)
    stop = timeit.default_timer()
    if LOG:
        logs_time.write("#registry_rnd_car#"+str(stop-start)+"\n")
    if VERBOSE:
        print("LP             :", lp)
        print("Qv            x:", Qv.toString().split()[0])
        print("              y:", Qv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("gQv           x:", gQv.toString().split()[0])
        print("              y:", gQv.toString().split()[1])


    return Qv,gQv


def com_setup_st_1_OBU(Qv, gQv, Qr, ec):
    s = random.getrandbits(256)
    t = random.getrandbits(256)
    y = random.getrandbits(256)
    
    
    yP = ec.mulJ(P, y)
    
    sgQv = ec.mulJ(gQv, s)
    Qv_byte = point_to_byte(Qv)
    sgQv_byte = point_to_byte(sgQv)
    yP_byte = point_to_byte(yP)

    
    start = timeit.default_timer()
    try:
        RVC.index(Qr.toString())
        RL = True
    except ValueError:
        RL = False
    if not RL:
        
        A1 = bilinear_map(gQv, Qr, ec)
        
        msg_byte = A1.real.to_bytes(64,"big")+t.to_bytes(64,"big")+Qv_byte+sgQv_byte+yP_byte

        M1 = BonehFranklin_block.encryp(BonehFranklin_block, msg_byte, Qr, gP, P, ec)

        stop = timeit.default_timer()

    if LOG:
        logs_time.write("#com_setup_st_1_OBU#"+str(stop-start)+"\n")
    msg = A1.toString()+" "+str(t)+" "+Qv.toString()+" "+sgQv.toString()+" "+yP.toString()
    if VERBOSE:
        print("s           =",s)
        print("---------------------------------------------------------------------------------------------------------------")
        print("A1          =",A1.toString())
        print("---------------------------------------------------------------------------------------------------------------")
        print("sgammaQv      x:",sgQv.toString().split()[0])
        print("              y:",sgQv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("msg         =",msg)

    return s, t, y, M1


def com_setup_st_2_RSU(xi, gQr, M1):
    start = timeit.default_timer()
    
    decM1 = BonehFranklin_block.decrypt(BonehFranklin_block, M1, gQr, ec)
    A1_obu = int.from_bytes(decM1[0:64],"big")
    t = int.from_bytes(decM1[64:128],"big")
    Qv = byte_to_point(decM1[128:256],ec)
    sgQv = byte_to_point(decM1[256:384],ec)
    yP = byte_to_point(decM1[384:],ec)

    validity = False
    
    try:
        RVC.index(Qv.toString())
        RL = True
    except ValueError:
        RL = False

    if not RL:
        A1_rsu = bilinear_map(Qv, gQr, ec)


        validity = A1_rsu.real == A1_obu and not RL

        xisgQv = ec.mulJ(sgQv, xi)
        
        
        xiQv = ec.mulJ(Qv, xi)
        to_aul  = bilinear_map(xiQv, yP, ec)
        txiQv = ec.mulJ(xiQv, t)
        #!!!!!!!!!!!!!! TA PART !!!!!!!!!!!!!!!!!!
        
        AUL.append(to_aul.HCpow(to_aul, beta, ec.q).real)
        
    stop = timeit.default_timer()
    if LOG:
        logs_time.write("#com_setup_st_2_RSU#"+str(stop-start)+"\n")
    if VERBOSE:
        print("Dec M1      = ", Qv.toString() + "/" + str(A1_obu) + "/"  + sgQv.toString())
        print("---------------------------------------------------------------------------------------------------------------")
        print("Rec Qv      = x:", Qv.toString().split()[0])
        print("              y:", Qv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("Rec A1      = ", A1_obu)
        print("")
        print("On revoc    =", RL)
        print("RSU A1      = ", A1_rsu.real)
        print("Validity    =", validity)
        print("---------------------------------------------------------------------------------------------------------------")
        print("xisgQv      = x:", xisgQv.toString().split()[0])
        print("              y:", xisgQv.toString().split()[1])
    if validity:
        return xisgQv, txiQv
    return False
    
def com_setup_st_3_OBU(xisgQv, txiQv, Qv, gQv, t, s):
    s_inv = Math.modular_inverse(s, ec.r)
    t_inv = Math.modular_inverse(t, ec.r)
    
    start = timeit.default_timer()
    
    xigQv = ec.mulJ(xisgQv, s_inv)
    xiQv = ec.mulJ(txiQv, t_inv)
    
    s1 = ec.add2(xigQv, xiQv)
    s2 = ec.add2(gQv, Qv)
     
    validity = bilinear_map(s1, P, ec).real == bilinear_map(s2, xi_P, ec).real
    stop = timeit.default_timer()
    if LOG:
        logs_time.write("#com_setup_st_3_OBU#"+str(stop-start)+"\n")
    if VERBOSE:
        print("xigQv       = x:", xigQv.toString().split()[0])
        print("              y:", xigQv.toString().split()[1])
        print("")
        print("xiQv       = x:", xiQv.toString().split()[0])
        print("              y:", xiQv.toString().split()[1])
        print("Validity:", validity)
    return xigQv, xiQv

def inc_report(xigQv, xiQv, y, M, T):
    
    
    a = random.getrandbits(256)
    a_inv = Math.modular_inverse(a, ec.r)
    
    A_ID = ec.mulJ(xiQv, a)
    
    A_1 = ec.mulJ(bP, a_inv)
    A_1 = ec.mulJ(A_1, y)
    
    A_TMP = ec.mulJ(xigQv, a)
    
    A_2 = ec.mulJ(gbP, a_inv)
    A_2 = ec.mulJ(A_2, y)
    
    
    start = timeit.default_timer()
    H = hashlib.sha256()
    H.update(A_ID.toString().encode())
    H.update(A_1.toString().encode())
    H.update(str(M).encode())
    H.update(str(T).encode())
    h = int(H.hexdigest(), 16)
    
    
    A_TMP = ec.mulJ(A_TMP, h)
    
    
    A_2 = ec.add2(A_TMP, A_2)
    
    
    
    stop = timeit.default_timer()
    if LOG:
        logs_time.write("#inc_report#"+str(stop-start)+"\n")
    if VERBOSE:
        print("A_ID        = x:", A_ID.toString().split()[0])
        print("              y:", A_ID.toString().split()[1])
        print("")
        print("A_1         = x:", A_1.toString().split()[0])
        print("              y:", A_1.toString().split()[1])
        print("")
        print("A_2         = x:", A_2.toString().split()[0])
        print("              y:", A_2.toString().split()[1])
        print("")
        print("M:", M)
        print("T:", T)
        
    return A_ID, A_1, A_2, M, T


def verify_inc_report(A_ID, A_1, A_2, M, T):
    start = timeit.default_timer()
    
    H = hashlib.sha256()
    H.update(A_ID.toString().encode())
    H.update(A_1.toString().encode())
    H.update(str(M).encode())
    H.update(str(T).encode())
    h = int(H.hexdigest(), 16)
    
    _A2 = ec.mulJ(A_ID, h)
    _A2 = ec.add2(_A2, A_1)
    
    
    MSG_INTEGRITY = bilinear_map(_A2, gP, ec).real == bilinear_map(A_2, P, ec).real
    
    try:
        AUL.index(int(bilinear_map(A_ID, A_1, ec).real))
        ELIGIBLE = True
    except ValueError:
        ELIGIBLE = False
    
    
    stop = timeit.default_timer()
    if LOG:
        logs_time.write("#verify_inc_report#"+str(stop-start)+"\n")
    if VERBOSE:
        print("A_ID        = x:", A_ID.toString().split()[0])
        print("              y:", A_ID.toString().split()[1])
        print("")
        print("A_1         = x:", A_1.toString().split()[0])
        print("              y:", A_1.toString().split()[1])
        print("")
        print("A_2         = x:", A_2.toString().split()[0])
        print("              y:", A_2.toString().split()[1])
        print("")
        print("M:", M)
        print("T:", T)
        print("INTEGRITY:", MSG_INTEGRITY)
        print("ELIGIBLE:", ELIGIBLE)
        
        
    return ELIGIBLE and MSG_INTEGRITY  
for _ in range(100000):
    ###############################################################
    ###################### TA Initialization ######################
    ###############################################################

    gamma = random.getrandbits(256)
    beta = random.getrandbits(256)


    ec = EC.EC()
    ec.basepoint = ec.at(None)

    P = point_from_string("5043129300441609945407506669182625058536335471322304775938584243304718261194640451519661940563248010695774457408090415592610241477449099234529292056928008 666247954662654463404525465677532364696032279767598862589447399628719603485134381749222527960645298941726840184644466632693772818656650944574674921112187")

    gP = ec.mulJ(P, gamma) # gammaP
    bP = ec.mulJ(P, beta) # betaP
    gbP = ec.mulJ(bP, gamma) # gamma beta P

    AUL = [] # Anonimized User List
    RVC = [] # Revocation List


    ###############################################################
    ################## Participant registrations ##################
    ###############################################################


    print("###############################################################")
    print("################## Participant registrations ##################")
    print("###############################################################")
    Qv,gQv = registry_rnd_car("IDX-411" + str(time.time()), gamma)
    print("###############################################################" )
    Qw,gQw = registry_rnd_car("IDX-412" + str(time.time()), gamma)
    print("###############################################################")
    Qr, gQr = registry_rnd_car("RSU" + str(time.time()), gamma)

    xi = random.getrandbits(256)
    xi_P = ec.mulJ(P, xi)

    ###############################################################
    ##################### Communication Setup #####################
    ###############################################################


    print("###############################################################")
    print("##################### Communication Setup #####################")
    print("###############################################################")
    print()
    print("########################## OBU step ###########################")
    print()

    s, t, y, M1 = com_setup_st_1_OBU(Qv, gQv, Qr, ec)


    print()
    print("########################## RSU step ###########################")
    print()

    xisgQv, txiQv = com_setup_st_2_RSU(xi, gQr, M1)


    print()
    print("########################## OBU step ###########################")
    print()

    xigQv, xiQv = com_setup_st_3_OBU(xisgQv, txiQv, Qv, gQv, t, s)

    ###############################################################
    ####################### Incident Report #######################
    ###############################################################


    print("###############################################################")
    print("####################### Incident Report #######################")
    print("###############################################################")
    print()
    print("########################### Create ############################")
    print()


    A_ID, A_1, A_2, M, T = inc_report(xigQv, xiQv, y, "baleset", time.time())

    print()
    print("########################### Verify ############################")
    print()
    verify_inc_report(A_ID, A_1, A_2, M, T)