"""Render the sourced governance case timeline and its browser evidence table."""
from pathlib import Path
import json
import copy
import html
import re
from datetime import date
import calendar
import sys

UNGROUPED = False

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output'
OUT.mkdir(exist_ok=True)
data = json.loads((ROOT / 'sources/governance-combined-data.json').read_text())
evidence_events={lane['name']:copy.deepcopy(lane['events']) for lane in data['lanes']}
oversight={
    'Human-participant protections':[(1962.78,'binding'),(1964.46,'guidance')],
    'Biological weapons / BWC':[(1975.23,'binding')],
    'Pathogen research and controls':[(1974,'guidance'),(1997.285,'binding')],
    'Life sciences / journal safeguards':[(2003.13,'guidance')],
    'IVF / UK fertility rules':[(1985,'guidance'),(1991.58,'binding')],
    'Laboratory biosafety':[(1975.25,'binding'),(1983,'guidance')],
    'Gene synthesis / screening':[(2009.88,'guidance'),(2025.315,'binding')],
    'Engineered organisms':[(1975,'guidance'),(1976.48,'binding')],
    'Gene therapy':[(1990,'binding')],
    'Modified organisms / Cartagena':[(2003.69,'binding')],
    'Stem-cell research':[(1991.58,'binding'),(2006,'guidance')],
    'Engineered gene drives / guidance':[(2014.75,'binding'),(2016,'guidance')],
    'AI laboratories / biosecurity':[(2024.47,'guidance')],
    'CRISPR / heritable editing':[(2015.92,'guidance'),(2016.085,'binding')],
    'Human embryo models / code':[(2021.4,'guidance')],
}
first_oversight={name:min(year for year,kind in entries) for name,entries in oversight.items()}
first_oversight['Human-participant protections']=1947
# Sorting includes earlier oversight retained in the evidence, not only bar anchors.
ordering_audit = json.loads((ROOT / 'sources/governance-audited-data.json').read_text())
for lane in ordering_audit['lanes']:
    name = lane['legacy_name']
    dates = [float(record['valid_from'][:4]) for record in lane['governance'] + lane['restrictions']]
    if dates:
        first_oversight[name] = min(first_oversight.get(name, float('inf')), *dates)
for lane in data['lanes']:
    for event in lane.get('prior_governance', []):
        first_oversight[lane['name']] = min(first_oversight.get(lane['name'], float('inf')), event['year'])

ordering_notes = json.loads((ROOT / 'sources/governance-ordering.json').read_text())
for name, anchor in ordering_notes['additional_anchors'].items():
    first_oversight[name] = min(first_oversight.get(name, float('inf')), float(anchor['date'][:4]))
narrative_groups = [
    ('Research practice and safety', ['Human-participant protections', 'Laboratory biosafety', 'Life sciences / journal safeguards', 'AI laboratories / biosecurity']),
    ('Pathogens and security', ['Pathogen research and controls', 'Biological weapons / BWC']),
    ('Genetic engineering', ['Engineered organisms', 'Modified organisms / Cartagena', 'Gene synthesis / screening', 'Engineered gene drives / guidance', 'Gene therapy']),
    ('Embryos and reproduction', ['IVF / UK fertility rules', 'Stem-cell research', 'CRISPR / heritable editing', 'Human embryo models / code']),
]
# Include inherited legal frameworks and earlier evidence when ordering groups.
narrative_groups = [
    (title, sorted(names, key=lambda name: first_oversight.get(name, float('inf'))))
    for title, names in narrative_groups
]
lane_lookup = {lane['name']: lane for lane in data['lanes']}
lanes = [lane_lookup[name] for _, names in narrative_groups for name in names]
group_starts = {names[0]: title for title, names in narrative_groups}
if UNGROUPED:
    lanes.sort(key=lambda lane: (first_oversight.get(lane['name'], float('inf')), lane['name']))
    editing_lane = next(lane for lane in lanes if lane['name'] == 'CRISPR / heritable editing')
    lanes.remove(editing_lane)
    stem_index = next(index for index, lane in enumerate(lanes) if lane['name'] == 'Stem-cell research')
    lanes.insert(stem_index + 1, editing_lane)
    group_starts = {}
output_stem = 'governance-swimlanes-ungrouped' if UNGROUPED else 'governance-swimlanes'
shared_frameworks = [] if UNGROUPED else [
    dict(year=1975.23, lanes=['Biological weapons / BWC', 'Pathogen research and controls'], label='BWC: weapons prohibition', source='https://treaties.un.org/Pages/showDetails.aspx?objid=0800000280101653'),
    dict(year=1991.58, lanes=['IVF / UK fertility rules', 'Stem-cell research', 'CRISPR / heritable editing'], label='Embryo licensing [UK]', source='https://www.legislation.gov.uk/ukpga/1990/37/contents'),
    dict(year=1993.085, lanes=['Modified organisms / Cartagena', 'Engineered gene drives / guidance'], label='Shared GMO rules', source='https://publications.parliament.uk/pa/cm199899/cmhansrd/vo990310/text/90310w09.htm'),
]
shared_starts = {}
for framework in shared_frameworks:
    for name in framework['lanes']:
        shared_starts[name] = min(shared_starts.get(name, float('inf')), framework['year'])

split_data = json.loads((ROOT / 'sources/governance-split-bands.json').read_text())
splits = split_data['splits']
audited = json.loads((ROOT / 'sources/governance-audited-data.json').read_text())
collective_restrictions = {lane['legacy_name']:lane['restrictions'] for lane in audited['lanes']}
visual_review=json.loads((ROOT / 'sources/governance-visual-review.json').read_text())
audited['sources'].update(visual_review['sources'])
for name,records in visual_review['additional_restrictions'].items():
    collective_restrictions.setdefault(name,[]).extend(records)
meta = data['metadata']
meta['ordering']='Lanes are grouped by narrative theme: research practice and safety; pathogens and security; genetic engineering; embryos and reproduction. Within each group, lanes follow the earliest mapped oversight anchor, including earlier frameworks retained in the audit, nonbinding guidance and funding conditions, rather than only the displayed markers. These are overlapping topics, not exclusive scientific or legal categories.'
start_year = 1960
review_day=date.fromisoformat(split_data['review_cutoff'])
end_year=review_day.year+(review_day-date(review_day.year,1,1)).days/(366 if calendar.isleap(review_day.year) else 365)
# Keep earlier milestones in a separate note rather than plotting them outside the window.
early_events=[(lane['name'],event) for lane in lanes for event in lane['events'] if event['year'] < start_year]
for lane in lanes:
    early_events.extend((lane['name'],dict(e,label=e['label'])) for e in lane.get('prior_governance',[]))
lanes=[lane for lane in lanes if lane['name'] != 'Animal research / 3Rs and ASPA']
for lane in lanes:
    if lane['name']=='Modified organisms / Cartagena':
        lane['events']=[e for e in lane['events'] if e['label']!='2000 adopted']
    lane['events']=[event for event in lane['events'] if event.get('display',True) and event['year'] >= start_year and not (lane['name']=='Biological weapons / BWC' and 'treaty opened' in event['label'])]
shared_duplicates = {
    'Biological weapons / BWC': {'1975 entered into force'},
    'IVF / UK fertility rules': {'1991 HFEA licensing begins'},
    'Stem-cell research': {'1991 UK embryo licensing'},
}
for lane in lanes:
    lane['events'] = [event for event in lane['events'] if event['label'] not in shared_duplicates.get(lane['name'], set())]
for framework in shared_frameworks:
    lane=next(lane for lane in lanes if lane['name']==framework['lanes'][0])
    lane['events'].append(dict(year=framework['year'],kind='policy',label=framework['label'],overview_label=framework['label'],source=framework['source'],shared=True,detail='Shared framework across the connected lanes; see scope notes.'))
    lane['events'].sort(key=lambda event:event['year'])
shown={(lane['name'],event['label']) for lane in lanes for event in lane['events']}
for name,events in evidence_events.items():
    for event in events:
        event['selection_note']='Shown on chart' if (name,event['label']) in shown else 'Evidence only / earlier context'
export_data=dict(data,lanes=lanes,split_bands=split_data,collective_restrictions=collective_restrictions,evidence_events=evidence_events,visual_review=visual_review,ordering_notes=ordering_notes)
(OUT / f'{output_stem}-data.json').write_text(json.dumps(export_data, indent=2))

BG = '#ffffff'
INK = '#273039'
MUTED = '#66717a'
GRID = '#e9edef'
BAND = '#f0f2f3'
BAND_HEIGHT = 14
RED = '#c62828'
GREEN = '#197650'
BLUE = '#126d9b'
GOLD = '#8b713b'
EXPANSION = '#16877e'
OVERSIGHT = '#487ca3'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'path'})
fig, ax = plt.subplots(figsize=(22, 14), facecolor=BG)
ax.set_facecolor(BG)
fig.subplots_adjust(left=.145, right=.985, top=.840, bottom=.070)
ax.set_xlim(start_year-1, end_year+3)
# Allocate vertical space from measured labels, in physical points.
from matplotlib.font_manager import FontProperties
fig.canvas.draw()
probe=fig.canvas.get_renderer()
font=FontProperties(family='DejaVu Sans',size=12)
points_per_year=22*72*(.985-.145)/(end_year+3-(start_year-1))
positions=[]
band_heights=[]
lane_bottoms=[]
lane_tops=[]
event_rows={}
next_y=0
group_heading_positions=[]
for lane in lanes:
    if lane['name'] in group_starts:
        if group_heading_positions:
            next_y += 14
        group_heading_positions.append((group_starts[lane['name']],next_y+10))
        next_y += 26
    lane_tops.append(next_y)
    ends=[]
    earlier=[e for name,e in early_events if name==lane['name']]
    if earlier:
        early_text='; '.join(e['label'] for e in earlier)
        early_font=FontProperties(family='DejaVu Sans',size=7)
        width=probe.get_text_width_height_descent(early_text,early_font,False)[0]*72/probe.dpi
        ends=[float('-inf'),points_per_year+width+15]
    for index,event in enumerate(lane['events']):
        text=event.get('overview_label',event['label'])
        width=max(probe.get_text_width_height_descent(line,font,False)[0] for line in text.splitlines())*72/probe.dpi
        x=(event['year']-(start_year-1))*points_per_year
        left,right=x-width/2-5,x+width/2+5
        tier=next((j for j,end in enumerate(ends) if end<left),len(ends))
        if tier==len(ends):ends.append(right)
        else:ends[tier]=right
        event_rows[(lane['name'],index)]=tier
    tiers=len(ends)
    above=max(25,11+((max(0,tiers-1))//2)*17+14)
    below=max(25,11+((max(1,tiers-1)-1)//2)*17+14)
    center=next_y+above+4
    positions.append(center)
    band_heights.append(BAND_HEIGHT)
    lane_bottoms.append(center+below+4)
    next_y=lane_bottoms[-1]+10
header,footer=145,48
total_points=next_y+header+footer
fig.set_size_inches(22,total_points/72)
ax.set_position([.145,footer/total_points,.985-.145,next_y/total_points])
ax.set_ylim(next_y,0)
ax.set_xticks([1960,1970,1980,1990,2000,2010,2020,2026])
ax.xaxis.tick_top()
ax.tick_params(axis='x',length=0,pad=7,labelsize=10,labelcolor=MUTED)
ax.set_yticks([])
ax.grid(axis='x',color=GRID,linewidth=.6,zorder=0)
for spine in ax.spines.values():
    spine.set_visible(False)

# White background, red editorial accents, and neutral field lines.
fig.add_artist(Rectangle((.025,1-18/total_points),.035,4/total_points,transform=fig.transFigure,facecolor=RED,edgecolor='none'))
fig.text(.025,1-48/total_points,'The rules of discovery',fontsize=25,weight='bold',color=INK)
fig.text(.025,1-68/total_points,f'Scientific milestones, debates and governance responses · UK, US and international · 1960–2026',fontsize=11,color=INK)
handles=[Line2D([],[],marker='|',linestyle='none',color=MUTED,markersize=9,markeredgewidth=4.5,label='Event / rule / guidance'),
 Line2D([],[],marker='X',linestyle='none',color=RED,markersize=7,label='Adverse event / concern'),
 Line2D([],[],marker='D',linestyle='none',color=GREEN,markersize=6,label='Discovery / treatment'),
 Patch(facecolor='#f4d1d1',edgecolor='none',label='Restricted in at least one mapped setting'),
 Patch(facecolor='#d4e3ef',edgecolor='none',label='Binding oversight / funding conditions'),
 Patch(facecolor='#d5e8da',edgecolor='none',label='Nonbinding guidance / voluntary arrangements'),
 Patch(facecolor=BAND,edgecolor='none',label='Before mapped field-specific oversight')]
fig.legend(handles=handles[:3],loc='upper left',bbox_to_anchor=(.021,1-82/total_points),ncol=4,
           frameon=False,fontsize=9.5,handletextpad=.55,columnspacing=2.0)
fig.legend(handles=handles[3:],loc='upper left',bbox_to_anchor=(.021,1-105/total_points),ncol=4,
           frameon=False,fontsize=9.2,handletextpad=.55,columnspacing=2.0)
fig.text(.025,1-132/total_points,'FIELD / GROUP',ha='left',fontsize=8.5,color=MUTED,weight='bold')
for title, heading_y in group_heading_positions:
    y = fig.transFigure.inverted().transform(ax.transData.transform((start_year, heading_y)))[1]
    fig.text(.025, y, title, ha='left', va='center', fontsize=10, color=MUTED)
# Fine separators make each field easier to follow without adding boxes.
for index in range(len(positions)-1):
    boundary=lane_bottoms[index]+5
    ax.axhline(boundary,color='#f0f2f3',lw=.55,zorder=0)



short_names={
'Biological weapons / BWC':'Biological weapons',
'Human-participant protections':'Human research',
'Pathogen research and controls':'Pathogen research',
'Gene synthesis / screening':'Gene synthesis',
'Human research / Common Rule':'Human research · Common Rule',
'Animal research / 3Rs and ASPA':'Animal research · 3Rs / UK law',
'Life sciences / dual-use oversight':'Dual-use research · US oversight',
'Life sciences / journal safeguards':'Research publishing',
'Life sciences / China biosecurity law':'Biosecurity · China law',
'Modified organisms / EU release rules':'Modified organisms · EU rules',
'Modified organisms / Japan rules':'Modified organisms · Japan rules',
'Engineered gene drives / guidance':'Gene drives',
'Synthetic biology / US commission':'Synthetic biology · US commission',
'IVF / UK fertility rules':'IVF',
'Stem-cell research / embryos and iPSCs':'Stem cells · embryos / iPSCs',
'Human embryo models / code':'Embryo models',
'Pathogen research / US pause':'Pathogens · US funding pause',
'Pathogens / select-agent rules':'Pathogens · US select agents',
'Modified organisms / Cartagena':'Modified organisms',
'CRISPR / heritable editing':'Heritable editing',
'AI laboratories / biosecurity':'AI laboratories',
'Stem-cell research':'Stem-cell research',
}
replacements={
'Helsinki adoption':'Helsinki','entered into force':'in force',
'insulin approval':'insulin','Gelsinger death':'Gelsinger death',
'human lines':'hESC lines','voluntary pause':'pause',
'Cartagena Act effective':'Act effective','phased-testing guidance':'testing guidance',
'institutional review':'institutional review','publication safeguards':'editor safeguards',
'funding limit':'funding limit','rules replaced':'rules replaced',
'release directive':'release directive','commission report':'commission advice'
}
def exact_year(value, end=False):
    if len(value) == 4:
        return int(value) + (1 if end else 0)
    if len(value) == 7:
        year, month = map(int, value.split('-'))
        day = calendar.monthrange(year, month)[1] if end else 1
        point = date(year, month, day)
    else:
        point = date.fromisoformat(value)
    return point.year + (point-date(point.year,1,1)).days/(366 if calendar.isleap(point.year) else 365)


def scoped_span(record, center, height, restriction=False):
    left=max(start_year,exact_year(record['valid_from']))
    right=min(end_year,exact_year(record['valid_to'],True) if record['valid_to'] else exact_year(split_data['review_cutoff']))
    if right <= left:
        return
    color='#f4d1d1' if restriction else '#d5e8da' if record['authority_type']=='guidance' else '#d4e3ef'
    edge=RED if restriction else GREEN if record['authority_type']=='guidance' else BLUE
    supported=min(right,max(left,exact_year(record['supported_through'],True)))
    intervals=[(left,supported,False),(supported,right,True)] if record['status']!='documented' else [(left,right,False)]
    for x0,x1,uncertain in intervals:
        if x1<=x0: continue
        ax.add_patch(Rectangle((x0,center-height/2),x1-x0,height,
            facecolor='white' if uncertain else color,edgecolor=edge if uncertain or restriction else 'none',
            linewidth=.5 if uncertain else 0,linestyle=':' if uncertain else '-',hatch=None,zorder=4))
    if record['valid_to'] is None:
        ax.scatter([right],[center],marker='>',s=13,color=edge,zorder=5)

def mixed_band(split, center, height):
    """Combine named, scoped records without treating missing evidence as permission."""
    cutoff=min(end_year,exact_year(split_data['review_cutoff']))
    boundaries={start_year,cutoff}
    for row in split['rows']:
        for record in row['governance']+row['restrictions']:
            if record.get('display',True) is False: continue
            boundaries.add(max(start_year,min(cutoff,exact_year(record['valid_from']))))
            boundaries.add(max(start_year,min(cutoff,exact_year(record['valid_to'],True) if record['valid_to'] else cutoff)))
            boundaries.add(max(start_year,min(cutoff,exact_year(record['supported_through'],True))))
    points=sorted(boundaries)
    fills={'binding':'#d4e3ef','guidance':'#d5e8da','restricted':'#f4d1d1'}
    for left,right in zip(points,points[1:]):
        midpoint=(left+right)/2
        statuses=set(); uncertain=False
        for row in split['rows']:
            active=[]
            for restriction,records in [(False,row['governance']),(True,row['restrictions'])]:
                for record in records:
                    if record.get('display',True) is False: continue
                    end=exact_year(record['valid_to'],True) if record['valid_to'] else cutoff
                    if exact_year(record['valid_from']) <= midpoint < end:
                        status='restricted' if restriction else 'guidance' if record['authority_type']=='guidance' else 'binding'
                        active.append(status)
                        uncertain |= record['status']!='documented' and midpoint>exact_year(record['supported_through'],True)
            if 'restricted' in active: statuses.add('restricted')
            else: statuses.update(active)
        if not statuses: continue
        base='restricted' if 'restricted' in statuses else 'binding' if 'binding' in statuses else 'guidance'
        mixed=len(statuses)>1
        edge=RED if 'restricted' in statuses else BLUE
        ax.add_patch(Rectangle((left,center-height/2),right-left,height,
            facecolor=fills[base],edgecolor=edge,linewidth=0,hatch=None,zorder=4))
        if uncertain:
            ax.add_patch(Rectangle((left,center-height/2),right-left,height,
                facecolor='none',edgecolor=MUTED,linewidth=.65,linestyle=':',zorder=5))
    if any(r['valid_to'] is None for row in split['rows'] for r in row['governance']+row['restrictions'] if r.get('display',True)):
        ax.scatter([cutoff],[center],marker='>',s=15,color=MUTED,zorder=5)

annotations=[]
marker_records=[]
rendered_band_starts={}
for i,lane,band_height in zip(positions,lanes,band_heights):
    y=fig.transFigure.inverted().transform(ax.transData.transform((start_year,i)))[1]
    name=short_names.get(lane['name'],lane['name'].replace(' / ',' · '))
    split = splits.get(lane['name'])
    title_y = y
    fig.text(.025,title_y,name,ha='left',va='center',fontsize=12,weight='bold',color=INK)
    start=max(start_year,lane['start']) if lane['start'] is not None else start_year
    start=min(start,shared_starts.get(lane['name'],start))
    rendered_band_starts[lane['name']]=start
    ax.add_patch(Rectangle((start,i-band_height/2),end_year-start,band_height,facecolor=BAND,edgecolor='none',zorder=2))
    if lane['name'] in shared_starts:
        shared_start=max(start_year,shared_starts[lane['name']])
        ax.add_patch(Rectangle((shared_start,i-band_height/2),end_year-shared_start,band_height,
            facecolor='#d4e3ef',edgecolor='none',zorder=2.5))
    earlier=[e for name,e in early_events if name==lane['name']]
    early_names={'1947 Nuremberg':'Nuremberg 1947',
                 '1964 Helsinki adoption':'Helsinki 1964',
                 '1955 safety conferences':'Biosafety meetings · 1955',
                 '1969 human IVF':'Human IVF · 1969'}
    if earlier:
        text='; '.join(early_names.get(e['label'],e['label']) for e in earlier)
        history_label=ax.annotate('← '+text,(start_year,i),xytext=(5,-17),textcoords='offset points',
                    ha='left',va='center',fontsize=7,color=MUTED,zorder=8)
        if len(earlier)==1:
            history_label.set_url(earlier[0]['source'])
    elif False:
        ax.annotate(lane['origin_display'],(start,i),xytext=(-6,0),textcoords='offset points',
                    ha='right',va='center',fontsize=7,color=MUTED)
    if split:
        mixed_band(split,i,band_height)
    else:
        regimes=oversight.get(lane['name'],[])
        for year,kind in sorted(regimes,key=lambda entry:entry[1]=='binding'):
            ax.add_patch(Rectangle((year,i-band_height/2),end_year-year,band_height,
                facecolor='#d4e3ef' if kind=='binding' else '#d5e8da',edgecolor='none',zorder=3))
    # Red takes priority if any mapped setting restricts the specified activity.
    for record in collective_restrictions.get(lane['name'],[]):
        selected=visual_review['lane_scopes'].get(lane['name'],{}).get('overview_restrictions')
        if selected is not None and record['id'] not in selected:
            continue
        scoped_span(record,i,band_height,True)
    for event_index,event in enumerate(lane['events']):
        kind=event['kind']
        marker,color,size={'meeting':('s',MUTED,34),'success':('D',GREEN,48),'failure':('X',RED,72),'policy':('s',MUTED,34),'alternative':('D',BLUE,37),'milestone':('o',MUTED,34)}[kind]
        branch = None
        mark_y=i + ((-1 if branch==0 else 1)*band_height/4 if branch is not None else 0)
        if color == MUTED:
            if not event.get('shared'):
                ax.plot([event['year'],event['year']],[mark_y-4.5,mark_y+4.5],
                        color=MUTED,linewidth=4.5,solid_capstyle='butt',zorder=6)
        else:
            ax.scatter([event['year']],[mark_y],marker=marker,color=color,s=size,zorder=9 if kind=='failure' else 6,edgecolors=BG,linewidths=.8)
        label=event['label']
        for before,after in replacements.items(): label=label.replace(before,after)
        dy=-.27 if event.get('dy',-.25)<0 else .28
        display_label=re.sub(r'^\d{4}(?:[–-]\d{2,4})?\s+', '', label)
        # Keep labels centred on their year; use vertical placement for close events.
        if label=='2007 human iPSCs':
            dy=.8
        if 'UK research purposes expanded' in label:
            display_label='UK licences';dy=-.8
        if label=='2001 funding limit':
            display_label='US limit'
        if label=='2009 rules replaced':
            display_label='US revision'
        if label=='2016 testing guidance':
            dy=-.27
        if label=='1972 Tuskegee exposed':
            display_label='Tuskegee exposed';dy=-.60
        if label=='2018 births':
            dy=-.60
        if label=='2023 Casgevy*':
            dy=-.27
        label_names={
            'Nuremberg':'Nuremberg Code', 'Helsinki':'Helsinki',
            'Research Act':'US Research Act', 'Belmont':'Belmont principles',
            'Common Rule':'US Common Rule', 'UK trial rules':'UK trial regulations',
            'ASPA':'UK animal research law', 'treaty opened':'Opened for signature',
            'transfer controls':'US transfer controls', 'editor safeguards':'Journal safeguards',
            'federal review':'US federal review', 'institutional review':'Institutional review',
            'pause':'Voluntary pause', 'funding pause':'US funding pause',
            'P3CO review':'P3CO review framework', 'human IVF':'Human IVF',
            'IVF birth':'First IVF birth', 'Warnock':'Warnock report',
            'HFE Act':'UK fertility law', 'safety conferences':'Biosafety meetings',
            'WHO manual':'WHO biosafety manual', 'CDC/NIH BMBL':'US biosafety manual',
            'CDC anthrax incident':'CDC anthrax incident', 'WHO revision':'WHO manual updated',
            'synthetic gene':'First synthetic gene', 'IGSC protocol':'DNA screening protocol',
            'synthetic genome':'Synthetic genome', 'commission advice':'US ethics report',
            'NIH rules':'NIH DNA guidelines', 'insulin':'Insulin approved',
            'settlement':'Sponsor restrictions', 'hESC lines':'Human embryonic cells',
            'US limit':'US funding limits', 'US revision':'US funding expanded',
            'ISSCR guidelines':'ISSCR guidelines', 'proposal':'Gene-drive proposal',
            'testing guidance':'US testing guidance', 'WHO testing guidance':'WHO mosquito guidance',
            'Adam':'Adam robot', 'AIBL launch':'UK biosecurity lab',
            'births':'Edited births', 'patterning model':'Embryo model',
            'UK code':'UK embryo-model code',
        }
        display_label=label_names.get(display_label,display_label)
        if display_label=='revision':
            display_label='Common Rule updated' if lane['name']=='Human-participant protections' else 'ISSCR updated'
        if lane['name']=='AI laboratories / biosecurity' and display_label=='WHO guidance':
            display_label='WHO biosecurity guidance'
        if lane['name']=='CRISPR / heritable editing' and display_label=='summit':
            display_label='First summit' if event['year']<2018 else 'Second summit'
        if label=='1982 insulin': dy=-.60
        if label=='2001 funding limit': dy=.28
        if label=='2009 rules replaced': dy=.28
        if label=='2007 human iPSCs': dy=.85
        if label=='2006 mouse iPSCs': dy=-.27
        if label=='1998 hESC lines':
            display_label='Human ESCs';dy=-.27
        if label=='2009 rules replaced': display_label='US limits eased'
        if label=='2001 funding limit': display_label='US limits'
        if label=='2024 WHO guidance': display_label='WHO biosecurity'
        if label=='2026 AIBL launch': display_label='UK biosecurity lab'
        # Sentence case preserves established acronyms and proper names.
        if display_label and display_label[0].islower(): display_label=display_label[0].upper()+display_label[1:]
        rows=max(event_rows[(lane['name'],j)] for j in range(len(lane['events'])))+1
        row=event_rows[(lane['name'],event_index)]
        step=.44 if lane['name'] in {'Pathogen research and controls','Stem-cell research'} else .42
        dy=(-1 if row%2==0 else 1)*(.53+(row//2)*step)
        # Every event stays centred on the single band.
        if lane['name']=='Biological weapons / BWC':
            display_label='BWC in force'
        if lane['name']=='Modified organisms / Cartagena':
            display_label='Recombinant DNA' if event['kind']=='milestone' else 'Cartagena adopted' if event['year']<2001 else 'Cartagena in force' if event['year']<2004 else 'UK becomes party'
        if lane['name']=='Stem-cell research':
            display_label={
                'UK licences':'Research purposes expanded [UK]',
                'US limits':'Funding limits [US]',
                'US limits eased':'Funding expanded [US]',
                'ISSCR guidelines':'ISSCR guidance [Intl]',
                'ISSCR updated':'ISSCR updated [Intl]',
            }.get(display_label,display_label)
        display_label={'Gene-therapy proposal':'Initial proposal','First gene-therapy trial':'First therapeutic trial','Cartagena in force':'Cartagena in force','Editors’ biosecurity statement':'Editors’ statement','Select-agent transfer rules effective':'Transfer rules [US]','Possession rules effective':'Possession rules [US]','Institutional DURC implemented':'Institutional review [US]'}.get(display_label,display_label)
        display_label={'National Academies recommendations':'Academies recommendations','WHO mosquito guidance':'WHO mosquito guidance','UK embryo-editing licence':'Embryo-editing licence [UK]'}.get(display_label,display_label)
        if lane['name']=='Life sciences / journal safeguards': display_label='Joint statement / PNAS policy'
        if lane['name']=='IVF / UK fertility rules':
            display_label={'Voluntary licensing':'Voluntary IVF licensing','HFEA licensing begins':'Statutory IVF licensing','UK fertility law':'HFE Act passed'}.get(display_label,display_label)
        display_label={'UK trial reform effective':'UK trial reform\neffective','High-risk policy released':'High-risk policy','NIH procurement conditions':'NIH procurement\nconditions'}.get(display_label,display_label)
        concise={
            'US drug-trial safeguards':'Drug-trial rules',
            'US PHS review policy':'PHS review', 'US Research Act':'Research Act',
            'Belmont principles':'Belmont', 'US Common Rule':'Common Rule',
            'UK trial regulations':'UK trial rules', 'Revisions published':'Rule revised',
            'Revised Common Rule compliance':'New rule applies',
            'UK trial reform\neffective':'UK reform',
            'NIH-funded DNA conditions':'NIH guidelines',
            'GB workplace duties':'Workplace duties',
            'WHO biosafety manual':'WHO manual', 'US biosafety manual':'BMBL',
            'COSHH bio-agent duties':'COSHH',
            'CDC potential anthrax exposure':'CDC anthrax incident',
            'WHO manual updated':'WHO update',
            'Voluntary IVF licensing':'Voluntary licensing',
            'Statutory IVF licensing':'HFEA licensing',
            'First therapeutic trial':'First trial', 'Wilson sponsor ban':'Wilson ban',
            'Luxturna approved [US]':'Luxturna', 'Casgevy authorised':'Casgevy',
            'UK embryo licensing':'Embryo licensing [UK]',
            'Research purposes expanded [UK]':'Purposes expanded [UK]',
            'Funding limits [US]':'Funding limits [US]',
            'Funding expanded [US]':'Limits eased [US]',
            'ISSCR guidance [Intl]':'ISSCR guidelines', 'ISSCR updated [Intl]':'ISSCR update',
            'P3CO review framework':'P3CO',
            'Institutional review [US]':'Institutional review',
            'US federal review':'Federal review', 'US funding pause':'Funding pause',
            'NIH DGOF restrictions':'NIH restrictions', 'High-risk policy':'High-risk policy',
            'Joint statement / PNAS policy':'Editors’ statement / PNAS',
            'Cartagena in force':'Cartagena', 'UK becomes party':'UK joins',
            'Chemical gene synthesis':'Gene synthesised',
            'IGSC screening protocol':'IGSC screening',
            'UK screening guidance':'UK screening',
            'NIH procurement\nconditions':'NIH conditions',
            'Gene-drive proposal':'Burt proposal',
            'Academies recommendations':'Academies report',
            'GB contained-use law':'Contained-use rules [GB]',
            'WHO mosquito guidance':'WHO guidance',
            'Programmable CRISPR':'CRISPR-Cas9',
            'Embryo-editing licence [UK]':'Research licence [UK]',
            'Edited births reported':'Edited births',
            'Human patterning model':'Patterning model',
            'UK embryo-model code':'UK code',
            'ISSCR model guidance':'ISSCR guidance', 'ISSCR model update':'ISSCR update',
            'Robot scientist':'Robot scientist', 'WHO biosecurity':'WHO guidance',
        }
        display_label=concise.get(display_label,display_label)
        if display_label=='Wilson ban': display_label='Wilson sponsor limit'
        display_label=event.get('overview_label',display_label)
        ann=ax.annotate(display_label,xy=(event['year'],mark_y),xytext=(event['year'],i+dy),ha='center',va='center',fontsize=12,color=INK,zorder=7)
        ann.set_url(event['source'])
        annotations.append((ann,i, bool(split)))
        ann._layout_tier=event_rows[(lane['name'],event_index)]
        marker_records.append(dict(lane=lane['name'],label=event['label'],display_label=display_label,year=event['year'],row=branch,marker_y=mark_y,expected_y=i+((-1 if branch==0 else 1)*band_height/4 if branch is not None else 0)))

# A solid grey bar extends the marker between related rows.
connector_labels=[]
centres = {lane['name']: center for lane, center in zip(lanes, positions)}
for record in shared_frameworks:
    x=record['year']
    ys=[centres[name] for name in record['lanes']]
    ax.plot([x,x],[min(ys)-4.5,max(ys)+4.5],color=MUTED,linewidth=4.5,
            solid_capstyle='butt',zorder=5)
if shared_frameworks:
    note_y=(centres['IVF / UK fertility rules']+centres['Stem-cell research'])/2
    shared_note=ax.annotate('Indicating when a ruling\nimpacts multiple fields',
        xy=(1991.58,note_y),xytext=(1987,note_y),
        ha='right',va='center',fontsize=12,fontstyle='italic',color=MUTED,
        bbox=dict(facecolor=BG,edgecolor='none',pad=2),
        arrowprops=dict(arrowstyle='->',color=MUTED,lw=.8,
                        connectionstyle='arc3,rad=.15',shrinkA=5,shrinkB=5),
        zorder=10)
    connector_labels.append(shared_note)
    ipsc_event=next(event for lane in lanes if lane['name']=='Stem-cell research' for event in lane['events'] if event['label']=='2007 human iPSCs')
    ipsc_note=ax.annotate('Reduce the need for hESCs\nin many applications',
        xy=(ipsc_event['year'],centres['Stem-cell research']-23),
        xytext=(2017,centres['Stem-cell research']-52),
        ha='center',va='center',fontsize=11,fontstyle='italic',color=MUTED,
        bbox=dict(facecolor=BG,edgecolor='none',pad=2),
        arrowprops=dict(arrowstyle='->',color=MUTED,lw=.8,
                        connectionstyle='arc3,rad=.1',shrinkA=5,shrinkB=5),zorder=10)
    ipsc_note.set_url('https://research.nih.gov/research-at-nih/induced-pluripotent-stem-cell-research')
    connector_labels.append(ipsc_note)


fig.text(.025,28/total_points,'Red applies only to the named activity: hESC eligibility is not iPSC research; reproductive editing is not CRISPR research. No worldwide ban is implied.',fontsize=9,color=MUTED)
fig.text(.025,14/total_points,'Open arrows mean no verified end; dotted outlines mark uncertain continuity. Named groups, activities and audit sources are detailed below.',fontsize=8,color=MUTED)
fig.canvas.draw()
# Use the reserved tier within each lane; never move a label into another lane.
renderer=fig.canvas.get_renderer()
for ann,lane_y,is_split in annotations:
    tier=ann._layout_tier
    above=tier%2==0
    offset=11+(tier//2)*17
    ann.set_anncoords('offset points')
    ann.set_position((0,offset if above else -offset))
    ann.set_verticalalignment('bottom' if above else 'top')
    if shared_frameworks:
        ann.set_bbox(dict(facecolor=BG, edgecolor='none', pad=.6))
fig.canvas.draw()
for ann,i,is_split in annotations:
    box=ann.get_window_extent(renderer)
    if box.x0<0 or box.x1>fig.bbox.width:
        print('Label outside canvas:',i,ann.get_text())
lane_bounds={center:(top,bottom) for center,top,bottom in zip(positions,lane_tops,lane_bottoms)}
for ann,center,_ in annotations:
    box=ann.get_window_extent(renderer)
    top,bottom=lane_bounds[center]
    top_px=ax.transData.transform((start_year,top))[1]
    bottom_px=ax.transData.transform((start_year,bottom))[1]
    assert bottom_px <= box.y0 and box.y1 <= top_px, 'Label intrudes into another lane: '+ann.get_text()
    assert 0 <= box.x0 and box.x1 <= fig.bbox.width, 'Label clipped: '+ann.get_text()
assert all(m['marker_y']==m['expected_y'] for m in marker_records), 'Marker is off its assigned row centre.'
boxes=[ann.get_window_extent(renderer) for ann,_,_ in annotations]
for label in connector_labels:
    box=matplotlib.text.Text.get_window_extent(label,renderer)
    assert not any(box.overlaps(other) for other in boxes), 'Note overlaps: '+label.get_text()+' / '+str([ann.get_text() for ann,_,_ in annotations if box.overlaps(ann.get_window_extent(renderer))])
    assert 0 <= box.x0 and box.x1 <= fig.bbox.width, 'Shared-framework label clipped.'
overlap_count=sum(boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i))
assert overlap_count==0, f'{overlap_count} event-label overlaps remain.'
(OUT/('governance-ungrouped-validation.json' if UNGROUPED else 'governance-split-validation.json')).write_text(json.dumps(dict(event_count=len(marker_records),lane_order=[lane['name'] for lane in lanes],ordering_years=first_oversight,band_heights=band_heights,band_starts=rendered_band_starts,band_end=end_year,split_lane_count=len(splits),event_label_overlaps=int(overlap_count),lane_label_intrusions=0,markers=marker_records),indent=2))
for ext in ['svg','png','pdf']:
    fig.savefig(OUT/f'{output_stem}.{ext}',dpi=170,facecolor=BG)

if UNGROUPED:
    sys.exit(0)

rows=[]
for lane in lanes:
    rows.append(f'<tr><th>{html.escape(lane["name"])}</th><td>Origin</td><td>{html.escape(lane["origin"])}. {html.escape(lane.get("case_detail",""))}</td><td><a href="{lane["origin_source"]}">Source</a></td></tr>')
    for event in evidence_events[lane['name']]:
        rows.append(f'<tr><td>{html.escape(lane["name"])}</td><td>{event["year"]:.2f} · {event["kind"]} · {html.escape(event.get('selection_note','Historical detail'))}</td><td>{html.escape(event["label"])}. {html.escape(event["detail"])}</td><td><a href="{event["source"]}">Source</a></td></tr>')
    for interval in lane['restrictions']:
        rows.append(f'<tr><td>{html.escape(lane["name"])}</td><td>{html.escape(interval["label"])}</td><td>{html.escape(interval["detail"])}</td><td><a href="{interval["source"]}">Source</a></td></tr>')
page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>The rules of discovery — swimlanes</title><style>body{margin:0;background:#ffffff;color:#252525;font:16px Arial,sans-serif}main{padding:18px 20px;max-width:1900px;margin:auto}img{width:100%;height:auto}a{color:#c62828}details{margin:20px 0;padding:16px 0;border-top:1px solid #e4e8ea}summary{cursor:pointer;font-weight:600}table{border-collapse:collapse;width:100%;font-size:14px}td,th{text-align:left;vertical-align:top;padding:10px;border-bottom:1px solid #cfdbdf}p{max-width:1100px;line-height:1.5}</style><main><nav style="font-size:13px;margin-bottom:12px"><a href="governance-swimlanes.svg">Full-size figure</a> &nbsp; <a href="governance-swimlanes.pdf">PDF</a> &nbsp; <a href="governance-swimlanes.png">PNG</a></nav><img src="governance-swimlanes.svg" alt="15 UK, US and international case lanes on a shared 1960 to 2026 year axis, arranged in four light narrative groups, with heritable editing directly below stem-cell research. Solid grey bars mark events, rules and guidance; red crosses mark adverse events; green diamonds mark discoveries and treatments; pale red segments mark scoped restrictions."><details><summary>Sources, scope, and event details</summary>'''
for key in ['colour_rule','ordering','origin_rule','restriction_rule','success_rule','failure_rule','coverage']:
    page+='<p>'+html.escape(meta[key])+'</p>'
# Inline SVG preserves the source hyperlinks on event labels.
svg=(OUT/'governance-swimlanes.svg').read_text()
svg=svg[svg.index('<svg'):]
svg=svg.replace('<svg ', '<svg id="governance-chart" role="img" aria-label="Scientific governance timeline in four narrative groups. Event labels link to sources." style="display:block;width:100%;height:auto" ',1)
page=re.sub(r'<img src="governance-swimlanes.svg"[^>]*>',lambda match:svg,page,count=1)
page=page.replace('</style>', '#governance-chart a{cursor:pointer}#governance-chart a:hover{opacity:.65}#governance-chart a:focus-visible{outline:2px solid #126d9b}</style>')
page+='<p>Grey means before the mapped field-specific oversight, not an absence of regulation. General law, consent duties and professional ethics may already apply. Blue and green show selected frameworks within their stated jurisdictions and activities. This timeline does not establish the first oversight of every field worldwide.</p>'
page+='<p>Overview selection: supporting events remain in this table. Grouped entries retain the final named step as their marker date; related earlier steps are listed separately. Close dates do not imply identical events or causation.</p>'
page+='<p>Extended grey bars show shared frameworks, not identical rules or the existence of later techniques at the framework date. Connected bands extend to the shared legal framework; scientific milestones keep their own dates. The BWC prohibits weapons purposes, not peaceful pathogen research. Embryo licensing and GMO rules apply only within their defined subjects and activities.</p>'
page+='<h3>Lane ordering</h3><p>'+html.escape(ordering_notes['rule'])+'</p><ul>'
for name, anchor in ordering_notes['additional_anchors'].items():
    page+='<li>'+html.escape(name+' — '+anchor['date']+': '+anchor['label']+'. '+anchor['detail'])+' <a href="'+html.escape(anchor['source'],quote=True)+'">Source</a>. '+html.escape(anchor['reading_depth'])+'</li>'
page+='</ul>'
page+='<h3>Collective colour rule and scope</h3><p>Red takes priority whenever at least one mapped setting restricts specified work or an actor. This includes funding conditions, voluntary pauses and sponsor-specific restrictions; it is not a claim of a worldwide or whole-field legal ban. Otherwise binding requirements take priority over guidance. Grey means no framework mapped, not unrestricted.</p>'
for name,split in splits.items():
    page+='<h4>'+html.escape(name)+'</h4>'
    for row in split['rows']:
        page+='<p><strong>'+html.escape(row['label'])+'</strong> '+html.escape(row.get('note',''))+'</p><ul>'
        for record in row['governance']+row['restrictions']:
            links=' · '.join('<a href="'+html.escape(split_data['sources'][sid]['url'],quote=True)+'">Source</a> ('+html.escape(split_data['sources'][sid]['review_depth'])+')' for sid in record['source_ids'])
            page+='<li>'+html.escape(record['short_label']+': '+record['subjects']+' '+record['continuation_basis'])+' '+links+'</li>'
        page+='</ul>'
page+='<h3>Scope corrections and incident responses</h3>'
for name,scope in visual_review['lane_scopes'].items():
    page+='<p><strong>'+html.escape(name)+'</strong>: '+html.escape(scope['explanation'])+'</p>'
for finding in visual_review['findings']:
    page+='<p>'+html.escape(finding['text'])+' '+ ' · '.join('<a href="'+html.escape(visual_review['sources'][sid]['url'],quote=True)+'">'+html.escape(visual_review['sources'][sid]['title'])+'</a> ('+html.escape(visual_review['sources'][sid]['review_depth'])+')' for sid in finding['source_ids'])+'</p>'
page+='<h3>Restrictions and their exact scope (including audit-only records)</h3>'
for name,records in collective_restrictions.items():
    for record in records:
        links=' · '.join('<a href="'+html.escape(audited['sources'][sid]['url'],quote=True)+'">Source</a>' for sid in record['source_ids'])
        page+='<p><strong>'+html.escape(name+' — '+record['short_label'])+'</strong><br>'+html.escape(record['valid_from']+' to '+(record['valid_to'] or 'open; no verified end')+'. '+record['subjects']+' '+record['continuation_basis'])+' '+links+'</p>'
page+='<p><a href="../sources/2026-11-05-correction-review.md">Full audit details</a></p>'
page+='<h3>Before 1960: foundations</h3><p>The chart begins in 1960. Earlier milestones remain important context:</p><ul>'+''.join(f'<li><a href="{html.escape(e["source"],quote=True)}">{html.escape(e["label"])}</a> — {html.escape(name)}</li>' for name,e in sorted(early_events,key=lambda item:item[1]['year']))+'</ul>'
page+='<p>Sources accessed 5 November 2026. New evidence was reviewed as official indexed excerpts, metadata, abstracts, or selected webpage passages. Full legal instruments and research papers were not reviewed in this pass. The iPSC sources establish reprogramming milestones, not universal clinical safety. Decimal positions are approximate where only a year is known. This design is Economist-inspired; it is not an Economist publication.</p><table><thead><tr><th>Case</th><th>Date / mark</th><th>Meaning and limits</th><th>Evidence</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table></details></main></html>'
(OUT/'governance-swimlanes.html').write_text(page)
print(f'Rendered {len(lanes)} lanes and {sum(len(l["events"]) for l in lanes)} event markers.')
