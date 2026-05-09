import numpy as np
from scipy.stats import norm
from math import *

# buyer = [demand,std.dev,transportationCost,ordering cost,holding cost,backorder cost,lead time]
# vendor = [production rate,setup cost,holding cost]
buyer1 = []
buyer2 = []
buyer3 = []
buyer4 = []
buyer5 = []
buyers = [buyer1,buyer2,buyer3,buyer4,buyer5]

def vendorcost(vendor,n,D,Q,sigma_di):
    part1 = (vendor[1]*D)/(n*Q)
    part2 = (2-n)*sigma_di*vendor[0]+n-1
    part3 = vendor[2]*(sigma_di/2)*(Q/D)

    total = part1 + part2*part3
    return total


def psi(ki):
    return norm.pdf(ki)-ki*(1-norm.cdf(ki))


def expected_shortage(buyer,ki):
    return buyer[1]*sqrt(buyer[6])*psi(ki)


def buyer_cost(buyer,ki,n,Qi):
    part1 = ((buyer[3]/n)+buyer[2])*(buyer[0]/Qi)
    part2 = buyer[4]*(Qi/2+ki*buyer[1]*sqrt(buyer[6]))
    part3 = buyer[5]*buyer[0]*expected_shortage(buyer,ki)/Qi

    return part1 + part2 + part3


def integrated_total_cost(vendor,buyers,n,i,Q,noOfbuyers,ki_values):
    sigma_di = 0
    Buyer = 0
    for buyer in buyers:
        sigma_di += buyer[0]

    vendor = vendorcost(vendor,n,buyers[i][0],Q,sigma_di)
    Qi_values = [noOfbuyers]
    Qi_values[i] = Q
    for a in range(0,noOfbuyers):
        if a != i:
            Qi = Q*buyers[a][0]/buyers[i][0]
            Qi_values[a] = Qi

    for buyer in buyers:
        buyercost = buyer_cost(buyer,ki_values[i],n,Qi_values[i])
        Buyer += buyercost

    return vendor + buyer

def Q_dssh(buyers,vendor,n,i):

    part1 = 0
    part2 = 0
    sigma_di = 0
    di = buyers[i][0]
    for buyer in buyers:
        part1 += buyer[3]/n + buyer[2]
        part2 += buyer[4]*buyer[0]
    for buyer in buyers:
        sigma_di += buyer[0]
    An = vendor[2]*((2-n)*sigma_di*vendor[0]+n-1)*sigma_di/2*1/di

    Qsquare = 2*di*(part1+vendor[1]/n)/(part2+2*di*An)

    return round(sqrt(Qsquare))

def Q(buyers,vendor,n,i,ki_values):
    part1 = 0
    part2 = 0
    part3 = 0
    sigma_di = 0
    b = 0
    di = buyers[i][0]
    for buyer in buyers:
        part1 += buyer[3] / n + buyer[2]
        part2 += buyer[4] * buyer[0]
        part3 += buyer[5]*expected_shortage(buyer,ki_values[b])
        b +=1
    for buyer in buyers:
        sigma_di += buyer[0]
    An = vendor[2] * ((2 - n) * sigma_di * vendor[0] + n - 1) * sigma_di / 2 * 1 / di

    Qsquare = 2 * di * (part1 + part3 + vendor[1] / n) / (part2 + 2 * di * An)

    return round(sqrt(Qsquare))

def get_ki(buyer,Q,D):
    p = 1-buyer[4]*Q/(buyer[5]*D)
    z_value = norm.ppf(p)

    return z_value

def kivalues(buyers,Q,D):
    values = []
    for buyer in buyers:
        values.append(get_ki(buyer,Q,D))

    return values

