from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
import struct
import os
import plotly.express as px
import plotly.graph_objs as go
import plotly.offline as pyo
import pandas as pd
print(os.getcwd())
os.chdir("/home/fduelli/customASO_project")
print(os.getcwd())
fileName = input("Please type the desired name for the graph and date files (no spaces): ")
targetExon = input("Please input the full sequence of the target exon to be changed for this ASO? ")
run = open("ASOgen.dat", "rb")
numAso_bytes = run.read(4)
print(numAso_bytes)
numAso = struct.unpack('>i', numAso_bytes)[0]
print(numAso)
length_bytes = run.read(4)
print(length_bytes)
length = struct.unpack('>i', length_bytes)[0]
print(length)
##start_index_bytes = run.read(4)
##start_index = struct.unpack('>i', start_index_bytes)[0]
import csv
with pd.ExcelWriter('ASO_walk_test_' + fileName + '.xlsx', engine='openpyxl') as writer:
    asoProg = 0
    allTargAcceptors = [] #Acceptor probabilities in the target exon
    allTargDonors = [] #Donor probabilities in the target exon
    allAsoArray = [] #Array saving each aso number increasing i.e ASO 1, ASO 2, ASO 3, ...
    allDonorPos = [] #Positions of each acceptor  in target exon  |
    allAcceptorPos = [] #Positions of each donor in target exon  -| can these two be combined?
    exonInd = 0 #starting index of the target exon
    untSeqAcceptors = [] #untreated sequence acceptor probabilities
    untSeqDonors = [] #untreated sequence donor probabilities
    for x in range(numAso+1):
        acceptor_normalized = [] #normalized acceptor probabilities based on difference from untreated sequence
        donor_normalized = [] #normalized donor probabilities based on difference from untreated sequence
        acceptor_colors = [] #graph colors for acceptor probabilities
        donor_colors = [] #graph colors for donor probabilities
        input_sequence = ''
        for _ in range(length):
            char = struct.unpack('>H', run.read(2))[0]
            input_sequence += chr(char)  # Convert int to char and add to the sequence string
        if asoProg == 0:
            exonInd = input_sequence.find(targetExon) - 10 #sets index of the target exon with the unedited gene sequence
        print(input_sequence)
        print(exonInd)
        context = 10000
        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        models = [load_model(resource_filename('spliceai', x)) for x in paths]
        x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
        y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

        acceptor_prob = y[0, :, 1]
        donor_prob = y[0, :, 2]

        #Print out the probabilities to console
        print("Acceptor probabilities:", acceptor_prob)
        print("Donor probabilities:", donor_prob)

        position = np.arange(0, len(acceptor_prob))
        if asoProg == 0: #saves the probabilities of acceptors and donors for the untreated sequence
            untSeqAcceptors = acceptor_prob
            untSeqDonors = donor_prob
        else: #saves the differences in probabilities between the current ASO and the untreated sequence
            tempCounter = 0
            for _ in untSeqAcceptors:
                acceptor_normalized.append(acceptor_prob[tempCounter] - untSeqAcceptors[tempCounter])
                donor_normalized.append(donor_prob[tempCounter] - untSeqDonors[tempCounter])
                tempCounter = tempCounter+1

        # save the probabilities to a CSV file
        data = {
        "position":position,
        "Acceptor": acceptor_prob,
        "Donor": donor_prob,
        }
        if(asoProg > 0):
            normData = {
            "position":position,
            "Normalized Acceptor": acceptor_normalized,
            "Normalized Donor": donor_normalized
            }
            normdf = pd.DataFrame(normData)
        # Create a DataFrame
        df = pd.DataFrame(data)
        # Write to an excel sheet
        if(asoProg == 0):
            df.to_excel(writer, sheet_name='no_ASO', index=False)
        else:
            df.to_excel(writer, sheet_name='Aso' + str(asoProg), index=False)
            normdf.to_excel(writer, sheet_name='ASO' + str(asoProg) + '(difference from untreated)', index=False)

        # Create color arrays with the same length as the probabilities, initialize with default color
        acceptor_colors.extend(['blue'] * len(acceptor_prob))
        donor_colors.extend(['orange'] * len(donor_prob))

        # Checks if the normalized values are greater than 10 percent of the total untreated probability, if yes checks if the acceptor or donor normalized values are positive or negative and displays green or red respectively
        tempCounter = 0
        if asoProg > 0:
            for _ in acceptor_normalized:
                if((abs(acceptor_normalized[tempCounter]) >= untSeqAcceptors[tempCounter]*0.1 and acceptor_normalized[tempCounter] >= 0.001) or (abs(donor_normalized[tempCounter]) >= untSeqDonors[tempCounter]*0.1 and donor_normalized[tempCounter] >= 0.001)):
                    if(acceptor_normalized[tempCounter] < 0):
                        acceptor_colors[tempCounter] = 'red'
                    elif(donor_normalized[tempCounter] < 0):
                        donor_colors[tempCounter] = 'red'     
                    elif(acceptor_normalized[tempCounter] > 0):
                        acceptor_colors[tempCounter] = 'green'
                    elif(donor_normalized[tempCounter] > 0):
                        donor_colors[tempCounter] = 'green'
                tempCounter = tempCounter+1

        # Create a line plot based on probabilities
        trace1 = go.Scatter(
            x=df['position'],
            y=df['Acceptor'],
            mode='lines+markers',
            name='Splice Acceptor Probability',
            line=dict(color='blue'),
            marker=dict(color=acceptor_colors)
        )
        trace2 = go.Scatter(
            x=df['position'],
            y=df['Donor'],
            mode='lines+markers',
            name='Splice Donor Probability',
            line=dict(color='orange'),
            marker=dict(color=donor_colors)
        )
        layout = go.Layout(
            title='Splice donor and acceptor probabilities for ' + str(asoProg),
            xaxis=dict(title='Position'),
            yaxis=dict(title='Probability')
        )
        fig = go.Figure(data=[trace1, trace2], layout=layout)
        pyo.plot(fig, filename='acceptor_donor_probabilities_'+ fileName + '_ASO' + str(asoProg) + '.html', auto_open=False)

        # Calculate target exon acceptor and donor probabilities and positions
        targExonAcceptor = acceptor_prob[exonInd: exonInd + len(targetExon) + 20].max()
        targExonDonor = donor_prob[exonInd: exonInd + len(targetExon) + 20].max()
        targExonAcceptorPos = exonInd + np.argmax(acceptor_prob[exonInd: exonInd + len(targetExon) + 20])
        targExonDonorPos = exonInd + np.argmax(donor_prob[exonInd: exonInd + len(targetExon) + 20])

        # Append probabilities and positions to lists
        allTargAcceptors.append(targExonAcceptor)
        allTargDonors.append(targExonDonor)
        allAsoArray.append(asoProg)
        allAcceptorPos.append(targExonAcceptorPos)
        allDonorPos.append(targExonDonorPos)

        # Print values for the greatest exon donors and acceptors in the target exon sequences
        print(targExonDonor)
        print(targExonAcceptor)
        asoProg += 1

    # Translate all arrays to numpy arrays for graphing
    allTargAcceptors = np.array(allTargAcceptors)
    allTargDonors = np.array(allTargDonors)
    allAsoArray = np.array(allAsoArray)
    allAcceptorPos = np.array(allAcceptorPos)
    allDonorPos = np.array(allDonorPos)

    # Compile results into a Dataframe
    allAsoData = {
        "ASOnum": allAsoArray,
        "Acceptor": allTargAcceptors,
        "Donor": allTargDonors,
        "Acceptor Position": allAcceptorPos,
        "Donor Position": allDonorPos
    }
allAsosdf = pd.DataFrame(allAsoData)
print(allAsosdf)
# Melt the DataFrame
long_Asosdf = allAsosdf.melt(id_vars=["ASOnum"], value_vars=["Acceptor", "Donor"], var_name="Type", value_name="Value")

# Add positions based on Type
long_Asosdf["Position"] = long_Asosdf.apply(
    lambda row: allAsosdf.loc[allAsosdf["ASOnum"] == row["ASOnum"], "Acceptor Position" if row["Type"] == "Acceptor" else "Donor Position"].values[0],
    axis=1,
)

# Separate DataFrames for Acceptor and Donor for creating specific traces
acceptor_df = long_Asosdf[long_Asosdf['Type'] == 'Acceptor']
donor_df = long_Asosdf[long_Asosdf['Type'] == 'Donor']

# Add the custom hover text column
acceptor_df["Hover Position"] = acceptor_df.apply(lambda row: f"{row['Type']} Position: {row['Position']}", axis=1)
donor_df["Hover Position"] = donor_df.apply(lambda row: f"{row['Type']} Position: {row['Position']}", axis=1)

# Create bar chart with separate traces
fig = go.Figure()

fig.add_trace(
    go.Bar(
        x=acceptor_df['ASOnum'], 
        y=acceptor_df['Value'], 
        name='Acceptor', 
        marker_color='blue',
        customdata=np.stack((acceptor_df['Type'], acceptor_df['Hover Position']), axis=-1),
        hovertemplate='<b>ASO number:</b> %{x}<br>' +
                      '<b>Type:</b> %{customdata[0]}<br>' +
                      '<b>Value:</b> %{y}<br>' +
                      '%{customdata[1]}<extra></extra>'
    )
)

fig.add_trace(
    go.Bar(
        x=donor_df['ASOnum'], 
        y=donor_df['Value'], 
        name='Donor', 
        marker_color='red',
        customdata=np.stack((donor_df['Type'], donor_df['Hover Position']), axis=-1),
        hovertemplate='<b>ASO number:</b> %{x}<br>' +
                      '<b>Type:</b> %{customdata[0]}<br>' +
                      '<b>Value:</b> %{y}<br>' +
                      '%{customdata[1]}<extra></extra>'
    )
)

# Update layout with title and axes labels
fig.update_layout(
    barmode='group',
    title="Acceptor and Donor probabilities for target Exon per ASO",
    xaxis_title="ASO number",
    yaxis_title="Acceptor and Donor Probabilities"
)

# Render the plot
pyo.plot(fig, filename="compiled_ASO_walk_results" + fileName + ".html", auto_open=False)

    

