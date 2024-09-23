import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics

target = "AKT1"
target_mdrun = "akt1_bbr"
rocs_run = "rocs_run2_0127"
nth = "2"
ROC_curve_name = "AKT1 bbr"

file_3cff= "/home/yeonji/Project/dude/" + target + "/" + rocs_run + "/" + target_mdrun + "_run" + nth + "_ROCS_RESULTS_1.rpt"
result_3cff = open(file_3cff, "r")
cont_3cff = result_3cff.readlines()
cont_3cff = cont_3cff[1:]

file_1cff = "/home/yeonji/Project/dude/" + target + "/" + rocs_run + "/" + target_mdrun + "_run" + nth + "_ROCS_1CFF_RESULTS_1.rpt"
result_1cff = open(file_1cff, "r")
cont_1cff = result_1cff.readlines()
cont_1cff = cont_1cff[1:]

cutoff_1cff = 0.003

perc_3cff = 0.01
perc_1cff = 0.01


##### Generate DataFrame #####

dic_3cff = {}
for line in cont_3cff:
    line_3cff = line.split()
    name_3cff = line_3cff[0]
    score_3cff = float(line_3cff[5])
    if name_3cff not in dic_3cff.keys():
        dic_3cff[name_3cff] = score_3cff
    elif name_3cff in dic_3cff.keys():
        score = dic_3cff[name_3cff]
        if score < score_3cff:
            dic_3cff[name_3cff] = score_3cff
        elif score >= score_3cff:
            dic_3cff[name_3cff] = score

dic_1cff = {}
for line in cont_1cff:
    line_1cff = line.split()
    name_1cff = line_1cff[0]
    score_1cff = float(line_1cff[5])
    if score_1cff > cutoff_1cff:
        if name_1cff not in dic_1cff.keys():
            dic_1cff[name_1cff] = score_1cff
        elif name_1cff in dic_1cff.keys():
            score = dic_1cff[name_1cff]
            if score < score_1cff:
                dic_3cff[name_1cff] = score_1cff
            elif score >= score_1cff:
                dic_1cff[name_1cff] = score

df_3cff = pd.DataFrame(columns=["Name", "ColorTanimoto"])
df_3cff["Name"] = list(dic_3cff.keys())
df_3cff["ColorTanimoto"] = list(dic_3cff.values())
df_1cff = pd.DataFrame(columns=["Name", "ColorTanimoto"])
df_1cff["Name"] = list(dic_1cff.keys())
df_1cff["ColorTanimoto"] = list(dic_1cff.values())

df_3cff = df_3cff.sort_values(by="ColorTanimoto", ascending=False)
df_3cff = df_3cff.reset_index(drop=True)
# print(df_3cff)
df_1cff = df_1cff.sort_values(by="ColorTanimoto", ascending=False)
df_1cff = df_1cff.reset_index(drop=True)
print(df_1cff[:20])

##### ROC-based Enrichment Factor #####

num_decoy_3cff = 0
num_active_3cff = 0
for name in df_3cff["Name"]:
    if "ZINC" in name:
        num_decoy_3cff+=1
    else:
        num_active_3cff+=1
one_perc_3cff = round(perc_3cff * num_decoy_3cff)

one_perc_list_3cff = []
count = 0
for name in df_3cff["Name"]:
    if count < one_perc_3cff:
        one_perc_list_3cff.append(name)
        if "ZINC" in name:
            count += 1
    elif count >= one_perc_3cff:
        break

a_3cff = 0
for lig in one_perc_list_3cff:
    if "ZINC" not in lig:
        a_3cff+=1
rocE_3cff = round((a_3cff/num_active_3cff)/perc_3cff, 3)

num_decoy_1cff = 0
num_active_1cff = 0
for name in df_1cff["Name"]:
    if "ZINC" in name:
        num_decoy_1cff+=1
    else:
        num_active_1cff+=1
one_perc_1cff = round(perc_1cff * num_decoy_1cff)

one_perc_list_1cff = []
count = 0
for name in df_1cff["Name"]:
    if count < one_perc_1cff:
        one_perc_list_1cff.append(name)
        if "ZINC" in name:
            count += 1
    elif count >= one_perc_1cff:
        break

a_1cff = 0
for lig in one_perc_list_1cff:
    if "ZINC" not in lig:
        a_1cff+=1
rocE_1cff = round((a_1cff/num_active_1cff)/perc_1cff, 3)

print(a_3cff, a_1cff)
print(rocE_3cff, rocE_1cff)

##### Plot ROC curve #####

def indicator(value):
    if "ZINC" in value:
        return 0
    else:
        return 1


df_3cff["Name"] = df_3cff.apply(lambda row:indicator(row["Name"]), axis=1)
y_3cff = df_3cff["Name"]
score_3cff = df_3cff["ColorTanimoto"]
fpr, tpr, thresholds = metrics.roc_curve(y_3cff, score_3cff, pos_label=1)
c = metrics.auc(fpr, tpr)

plt.figure()
lw = 1.5
plt.plot(fpr, tpr, color='blue',
         lw=lw, label='ROC curve (area = %0.2f)' % c)
plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(ROC_curve_name + " 3colorff")
plt.legend(loc="lower right")
plt.show()

df_1cff["Name"] = df_1cff.apply(lambda row:indicator(row["Name"]), axis=1)
y_1cff = df_1cff["Name"]
score_1cff = df_1cff["ColorTanimoto"]
fpr, tpr, thresholds = metrics.roc_curve(y_1cff, score_1cff, pos_label=1)
c = metrics.auc(fpr, tpr)

plt.figure()
lw = 1.5
plt.plot(fpr, tpr, color='blue',
         lw=lw, label='ROC curve (area = %0.2f)' % c)
plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(ROC_curve_name + " 1colorff")
plt.legend(loc="lower right")
plt.show()