<template>
    <div class="d-flex flex-column p-5">
        <LoadingOverlay name="plot" class="w-100 h-100">
            <img class="w-100" v-if="plotImageUrl != null" :src="plotImageUrl"/>
            <img class="w-100" v-else :src="defaultImage" style="opacity: 30%"/>
        </LoadingOverlay>
        <div class="d-flex flex-column">
            <b-form-group
                label="Select Plot Parameters"
                label-for="plotParamInput"
            >
                <b-form-select 
                    v-model="fileTypeToPlot"
                    size="sm" 
                    class="mt-3"
                >
                    <template #first>
                        <!-- <b-form-select-option :value="null" disabled>-- Select a file type to plot --</b-form-select-option> -->
                    </template>

                    <b-form-select-option v-for="plottable in plottableFilenames" :key="plottable" :value="plottable">
                        {{ plottable }}
                    </b-form-select-option>
                </b-form-select>
            </b-form-group>
            <div v-if="fileTypeToPlot != null">
                <b-form-group
                    label="Select Plot Parameters"
                    label-for="plotParamInput"
                >
                    <b-form-select 
                        v-model="plotParameter" 
                        :options="['new', ...plotParameters]" 
                        size="sm" 
                        class="mt-3"
                    ></b-form-select>
                </b-form-group>
                <div class="d-flex flex-row w-100">
                    <b-form-group label="x-axis" v-slot="{ ariaDescribedby }" class="mx-3" :disabled="plotParameter != 'new'">
                        <b-form-radio-group
                            id="radio-group-1"
                            v-model="xAxisIndex"
                            :options="xAxisOptions"
                            :aria-describedby="ariaDescribedby"
                            name="x-axis"
                            stacked
                        ></b-form-radio-group>
                    </b-form-group>
                    
                    <b-form-group label="y-axes" v-slot="{ ariaDescribedby }" class="mx-3" :disabled="plotParameter != 'new'">
                        <template v-slot:label>
                            Y_Axis
                            <b-button variant="link" v-b-modal.flags-modal size="sm">(help)</b-button>
                        </template>
                        <b-form-checkbox-group
                            id="y-axes-checkboxes"
                            v-model="yAxisIndeces"
                            :options="yAxisOptions"
                            :aria-describedby="ariaDescribedby"
                            name="y-axes"
                            stacked
                        ></b-form-checkbox-group>
                    </b-form-group>

                    <b-form-group v-slot="{ ariaDescribedby }" class="mx-3" :disabled="plotParameter != 'new'">
                        <template v-slot:label>
                            Flags
                            <b-button variant="link" v-b-modal.flags-modal size="sm">(help)</b-button>
                        </template>
                        <b-form-textarea
                            id="textarea"
                            v-model="flags"
                            placeholder="Enter flags here"
                            :aria-describedby="ariaDescribedby"
                            rows="6"
                            max-rows="6"
                        ></b-form-textarea>
                    </b-form-group>
                </div>
                <b-button variant="primary" @click="fetchPlot" :disabled="!canCreate">Create Plot</b-button>
                <div class="m-3"  v-for="[prefix, name] in plottablePrefixNames" v-bind:key="prefix">
                    <b>{{ prefix }}:</b> {{ name }}
                </div>
                <b-form-group 
                    v-if="fileTypeToPlot in plottableFiles && plottableFiles[fileTypeToPlot].length > 1" 
                    label="Plotted Files" 
                    v-slot="{ ariaDescribedby }" 
                    class="m-3"
                >
                    <b-form-checkbox-group
                        id="plotted-files-checkboxes"
                        v-model="plottedFiles"
                        :aria-describedby="ariaDescribedby"
                        name="plotted-files"
                        stacked
                    >
                        <b-form-checkbox 
                            v-for="plotFile in sortedPlotFiles" 
                            :key="plotFile" 
                            :value="plotFile"
                        >
                            {{ plotFile.name }}
                        </b-form-checkbox>
                    </b-form-checkbox-group>
                </b-form-group>
            </div>
            <b-modal id="flags-modal" size="lg" title="Graphing Flags" scrollable hide-footer>
                <b-table-simple>
                    <b-thead>
                        <b-tr>
                            <b-td>Flag</b-td>
                            <b-td>Description</b-td>
                        </b-tr>
                    </b-thead>
                    <b-tbody>
                        <!-- <b-tr>
                            <b-td>-label=PAR1,PAR2,...</b-td>
                            <b-td>parameter value(s) for PAR1,PAR2,...legend uses value(s) of linp name(s) containing PARi's</b-td>
                        </b-tr> -->
                        <b-tr>
                            <b-td>-showColumns</b-td>
                            <b-td>print the columns names (if found in file)</b-td>
                        </b-tr>
                        <!-- <b-tr>
                            <b-td>-peaks</b-td>
                            <b-td>mark photon peaks at n*omega-Ip-Up or (-2 Up for Helium)</b-td>
                        </b-tr> -->
                        <b-tr>
                            <b-td>-eV</b-td>
                            <b-td>transform x-axis as x -> x^2/2 *Rydberg</b-td>
                        </b-tr>
                        <!-- <b-tr>
                            <b-td>-polar</b-td>
                            <b-td>polar plot of 1d or 2d data</b-td>
                        </b-tr> -->
                        <b-tr>
                            <b-td>-compare</b-td>
                            <b-td>compare multiple 2d plots to first</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-linV</b-td>
                            <b-td>plot on linear value-axis (default is log)</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-linY</b-td>
                            <b-td>plot on linear y-axis (default is log)</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-logY</b-td>
                            <b-td>plot on logarithmic y-axis (default)</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-vrange=[vmin,vmax]</b-td>
                            <b-td>function value range for ALL plots</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-xrange=[xmin,xmax]</b-td>
                            <b-td>plot axis range for All plots</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-yrange=[ymin,ymax]</b-td>
                            <b-td>plot axis range for All plots</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-erange=[emin,emax]</b-td>
                            <b-td>range for plotting errors</b-td>
                        </b-tr>
                        <!-- <b-tr>
                            <b-td>-rmax=r</b-td>
                            <b-td>radius in polar plot</b-td>
                        </b-tr> -->
                        <b-tr>
                            <b-td>-equalAx</b-td>
                            <b-td>2d with commensurate axes</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-lineoutX=x{,w}</b-td>
                            <b-td>lineout of 2d plot nearest to x or sum [x-w/2,x+w/2], similar for Y</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-normalize{=x0}</b-td>
                            <b-td>scale to maximal value = 1 (nearest to x0)</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-symX12</b-td>
                            <b-td>symmetrize 2d plots by 1 &lt;--> 2 (special for He)</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-ratio</b-td>
                            <b-td>ratio of multiple plots to first</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-maxgraph</b-td>
                            <b-td>print location of graph's maximum</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-mingraph</b-td>
                            <b-td>print location of graph's minimum</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-sum</b-td>
                            <b-td>add files and columns single plot</b-td>
                        </b-tr>
                        <!-- <b-tr>
                            <b-td>-plotfile=file</b-td>
                            <b-td>specify plot file name</b-td>
                        </b-tr> -->
                        <b-tr>
                            <b-td>-batch</b-td>
                            <b-td>do keep window open</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-points</b-td>
                            <b-td>plot points rather than lines</b-td>
                        </b-tr>
                        <b-tr>
                            <b-td>-surface</b-td>
                            <b-td>do surface rather than contour plots</b-td>
                        </b-tr>
                        <!-- <b-tr>
                            <b-td>-figure=file.fig</b-td>
                            <b-td>read defaults and additional style options (will be superseeded by command line flags)</b-td>
                        </b-tr> -->
                    </b-tbody>
                </b-table-simple>
            </b-modal>
        </div>
    </div>
</template>

<script>
import LoadingOverlay from '../overlay/LoadingOverlay.vue'
import store from '@/store';
import { InputService, PlotService } from '@/service/epolyscat-service';
import { eventBus } from '@/event-bus';
import { plotObjects } from '@/fileData';

export default {
    components: { LoadingOverlay },
    props: ["selectedRuns"],
    store,
    data() {
        return {
            defaultImage: (
                "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAoAAAAHgCAYAAAA10dzkAAAAAXNSR0IArs4c6QAAAMJlWElmTU0AKgAAAAgABgESAAMAAAABAAEAAAEaAAUAAAABAAAAVgEbAAUAAAABAAAAXgEoAAMAAAABAAIAAAExAAIAAAAxAAAAZodpAAQAAAABAAAAmAAAAAAAAABkAAAAAQAAAGQAAAABTWF0cGxvdGxpYiB2ZXJzaW9uMy4zLjQsIGh0dHBzOi8vbWF0cGxvdGxpYi5vcmcvAAAAA6ABAAMAAAABAAEAAKACAAQAAAABAAACgKADAAQAAAABAAAB4AAAAAA7NYyVAAAACXBIWXMAAA9hAAAPYQGoP6dpAAAB62lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNi4wLjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyIKICAgICAgICAgICAgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIj4KICAgICAgICAgPHRpZmY6T3JpZW50YXRpb24+MTwvdGlmZjpPcmllbnRhdGlvbj4KICAgICAgICAgPHhtcDpDcmVhdG9yVG9vbD5NYXRwbG90bGliIHZlcnNpb24zLjMuNCwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy88L3htcDpDcmVhdG9yVG9vbD4KICAgICAgPC9yZGY6RGVzY3JpcHRpb24+CiAgIDwvcmRmOlJERj4KPC94OnhtcG1ldGE+ChshMWsAADQiSURBVHgB7d0L0FTlfT/wHxdBowErKASFGFsbNFptIVyczphWJpjYSWx0YhgTkTLatEpspGlAEf5J26GJtV7qhXGmGcdWKsWmtrGUDMU0sYV4QXPB29hOIggFJEbwEkBl//ucdNd94eVVzy7Zfd79nJllz55znnOe8/mdl/2+Z/ecd0ClOoSBAAECBAgQIECgawQGds2e2lECBAgQIECAAIFCQAB0IBAgQIAAAQIEukxAAOyygttdAgQIECBAgIAA6BggQIAAAQIECHSZgADYZQW3uwQIECBAgAABAdAxQIAAAQIECBDoMgEBsMsKbncJECBAgAABAgKgY4AAAQIECBAg0GUCAmCXFdzuEiBAgAABAgQEQMcAAQIECBAgQKDLBATALiu43SVAgAABAgQICICOAQIECBAgQIBAlwkIgF1WcLtLgAABAgQIEBAAHQMECBAgQIAAgS4TEAC7rOB2lwABAgQIECAgADoGCBAgQIAAAQJdJiAAdlnB7S4BAgQIECBAQAB0DBAgQIAAAQIEukxAAOyygttdAgQIECBAgIAA6BggQIAAAQIECHSZgADYZQW3uwQIECBAgAABAdAxQIAAAQIECBDoMgEBsMsKbncJECBAgAABAgKgY4AAAQIECBAg0GUCAmCXFdzuEiBAgAABAgQEQMcAAQIECBAgQKDLBATALiu43SVAgAABAgQICICOAQIECBAgQIBAlwkIgF1WcLtLgAABAgQIEBAAHQMECBAgQIAAgS4TEAC7rOB2lwABAgQIECAgADoGCBAgQIAAAQJdJiAAdlnB7S4BAgQIECBAQAB0DBAgQIAAAQIEukxAAOyygttdAgQIECBAgIAA6BggQIAAAQIECHSZgADYZQW3uwQIECBAgAABAdAxQIAAAQIECBDoMgEBsMsKbncJECBAgAABAgKgY4AAAQIECBAg0GUCAmCXFdzuEiBAgAABAgQEQMcAAQIECBAgQKDLBATALiu43SVAgAABAgQICICOAQIECBAgQIBAlwkIgF1WcLtLgAABAgQIEBAAHQMECBAgQIAAgS4TEAC7rOB2lwABAgQIECAgADoGCBAgQIAAAQJdJiAAdlnB7S4BAgQIECBAQAB0DBAgQIAAAQIEukxAAOyygttdAgQIECBAgIAA6BggQIAAAQIECHSZgADYZQW3uwQIECBAgAABAdAxQIAAAQIECBDoMgEBsMsKbncJECBAgAABAgKgY4AAAQIECBAg0GUCg7tsf1u6u3v37o3NmzfHu9/97hgwYEBL121lBAgQIECAwMERqFQq8dJLL8WYMWNi4MDuPBcmADZxbKXwN3bs2CbWoCkBAgQIECDQLoGNGzfGcccd167Nt3W7AmAT/OnMXxrSATRs2LAm1qQpAQIECBAg8IsS2LlzZ3ECp/Y+/ovabidtRwBsohq1j31T+BMAm4DUlAABAgQItEGg9j7ehk23fZPd+cF329l1gAABAgQIECDQPgEBsH32tkyAAAECBAgQaIuAANgWdhslQIAAAQIECLRPQABsn70tEyBAgAABAgTaIiAAtoXdRgkQIECAAAEC7RMQANtnb8sECBAgQIAAgbYICIBtYbdRAgQIECBAgED7BATA9tnbMgECBAgQIECgLQICYFvYbZQAAQIECBAg0D4BAbB99rZMgAABAgQIEGiLgADYFnYbJUCAAAECBAi0T0AAbJ+9LRMgQIAAAQIE2iIgALaF3UYJECBAgAABAu0TEADbZ2/LBAgQIECAAIG2CAiAbWG3UQIECBAgQIBA+wQEwPbZ2zIBAgQIECBAoC0CAmBb2G2UAAECBAgQINA+AQGwffa2TIAAAQIECBBoi4AA2BZ2GyVAgAABAgQItE9AAGyfvS0TIECAAAECBNoiIAC2hd1GCRAgQIAAAQLtExAA22dvywQIECBAgACBtggIgG1ht1ECBAgQIECAQPsEBMD22dsyAQIECBAgQKAtAgJgW9htlAABAgQIECDQPgEBsH32tkyAAAECBAgQaIuAANgWdhslQIAAAQIECLRPQABsn70tEyBAgAABAgTaIiAAtoXdRgkQIECAAAEC7RPIKgDecsstcfzxx8ehhx4akydPjoceeqhPueXLl8f48eOL5U899dRYsWLFAZf/7Gc/GwMGDIgbbrjhgMuYQYAAAQIECBDoDwLZBMBly5bFlVdeGYsWLYpHH300TjvttJg+fXps27at1zqsWbMmZsyYEbNnz47HHnsszj333OKxfv36/Zb/p3/6p/jud78bY8aM2W+eCQQIECBAgACB/iaQTQD8q7/6q7jkkkti1qxZcfLJJ8eSJUviXe96V3zta1/rtSY33nhjnH322fGFL3whTjrppPjTP/3T+I3f+I24+eabeyy/adOmmDNnTtx1111xyCGH9JjnBQECBAgQIECgPwpkEQD37NkT69ati2nTptVrMHDgwOL12rVr69MaR9L0xuXTvHTGsHH5vXv3xmc+85kiJH7gAx9obG6cAAECBAgQINBvBQbnsGfbt2+PN954I0aNGtWju+n1U0891WNa7cWWLVt6XT5Nrw1f+cpXYvDgwfG5z32uNqnP5927d0d61IadO3fWRj0TIECAAAECBLIRyOIM4MHQTGcU08fEd9xxR3Hxx9vZxuLFi2P48OH1x9ixY99OM8sQIECAAAECBDpKIIsAOHLkyBg0aFBs3bq1B156PXr06B7Tai/S9L6Wf+CBB4oLSMaNG1ecBUxnAp999tmYO3ducaVxbT2Nz/Pnz48dO3bUHxs3bmycbZwAAQIECBAgkIVAFgFwyJAhMWHChFi9enUdNX1/L72eOnVqfVrjSJreuHyat2rVqvry6bt/P/jBD+J73/te/ZGuAk4XjXzzm99sXFV9fOjQoTFs2LAej/pMIwQIECBAgACBTASy+A5gsky3gJk5c2ZMnDgxJk2aVNyv75VXXimuCk7zL7roojj22GMjfUybhiuuuCLOPPPMuO666+Kcc86Ju+++Ox555JG4/fbbi/kjRoyI9Ggc0lXA6czh+9///sbJxgkQIECAAAEC/UogmwB4wQUXxPPPPx8LFy6MdCHH6aefHitXrqxf6LFhw4ZIVwbXhjPOOCOWLl0aCxYsiKuuuipOPPHEuPfee+OUU06pLeKZAAECBAgQINCVAgMq1aEr97wFO52uAk4XhaTvBaaPhg0ECBAgQIBA5wt4/45485RZ59dLDwkQIECAAAECBFogIAC2ANEqCBAgQIAAAQI5CQiAOVVLXwkQIECAAAECLRAQAFuAaBUECBAgQIAAgZwEBMCcqqWvBAgQIECAAIEWCAiALUC0CgIECBAgQIBATgICYE7V0lcCBAgQIECAQAsEBMAWIFoFAQIECBAgQCAnAQEwp2rpKwECBAgQIECgBQICYAsQrYIAAQIECBAgkJOAAJhTtfSVAAECBAgQINACAQGwBYhWQYAAAQIECBDISUAAzKla+kqAAAECBAgQaIGAANgCRKsgQIAAAQIECOQkIADmVC19JUCAAAECBAi0QEAAbAGiVRAgQIAAAQIEchIQAHOqlr4SIECAAAECBFogIAC2ANEqCBAgQIAAAQI5CQiAOVVLXwkQIECAAAECLRAQAFuAaBUECBAgQIAAgZwEBMCcqqWvBAgQIECAAIEWCAiALUC0CgIECBAgQIBATgICYE7V0lcCBAgQIECAQAsEBMAWIFoFAQIECBAgQCAnAQEwp2rpKwECBAgQIECgBQICYAsQrYIAAQIECBAgkJOAAJhTtfSVAAECBAgQINACAQGwBYhWQYAAAQIECBDISUAAzKla+kqAAAECBAgQaIGAANgCRKsgQIAAAQIECOQkIADmVC19JUCAAAECBAi0QEAAbAGiVRAgQIAAAQIEchIQAHOqlr4SIECAAAECBFogIAC2ANEqCBAgQIAAAQI5CQiAOVVLXwkQIECAAAECLRAQAFuAaBUECBAgQIAAgZwEBMCcqqWvBAgQIECAAIEWCAiALUC0CgIECBAgQIBATgICYE7V0lcCBAgQIECAQAsEBMAWIFoFAQIECBAgQCAnAQEwp2rpKwECBAgQIECgBQICYAsQrYIAAQIECBAgkJOAAJhTtfSVAAECBAgQINACAQGwBYhWQYAAAQIECBDISUAAzKla+kqAAAECBAgQaIGAANgCRKsgQIAAAQIECOQkIADmVC19JUCAAAECBAi0QEAAbAGiVRAgQIAAAQIEchIQAHOqlr4SIECAAAECBFogIAC2ANEqCBAgQIAAAQI5CQiAOVVLXwkQIECAAAECLRAQAFuAaBUECBAgQIAAgZwEBMCcqqWvBAgQIECAAIEWCAiALUC0CgIECBAgQIBATgICYE7V0lcCBAgQIECAQAsEBMAWIFoFAQIECBAgQCAnAQEwp2rpKwECBAgQIECgBQICYAsQrYIAAQIECBAgkJOAAJhTtfSVAAECBAgQINACAQGwBYhWQYAAAQIECBDISUAAzKla+kqAAAECBAgQaIGAANgCRKsgQIAAAQIECOQkIADmVC19JUCAAAECBAi0QCCrAHjLLbfE8ccfH4ceemhMnjw5HnrooT4Jli9fHuPHjy+WP/XUU2PFihX15V977bX44he/GGn64YcfHmPGjImLLrooNm/eXF/GCAECBAgQIECgPwpkEwCXLVsWV155ZSxatCgeffTROO2002L69Omxbdu2XuuyZs2amDFjRsyePTsee+yxOPfcc4vH+vXri+VfffXVYj3XXHNN8fz1r389nn766fjYxz7W6/pMJECAAAECBAj0F4EBleqQw86kM34f/OAH4+abby66u3fv3hg7dmzMmTMn5s2bt98uXHDBBfHKK6/EfffdV583ZcqUOP3002PJkiX1aY0jDz/8cEyaNCmeffbZGDduXOOsXsd37twZw4cPjx07dsSwYcN6XcZEAgQIECBAoLMEvH9HZHEG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gACYbel0nAABAgQIECBQTkAALOemFQECBAgQIEAgWwEBMNvS6TgBAgQIECBAoJyAAFjOTSsCBAgQIECAQLYCAmC2pdNxAgQIECBAgEA5AQGwnJtWBAgQIECAAIFsBQTAbEun4wQIECBAgACBcgICYDk3rQgQIECAAAEC2QoIgNmWTscJECBAgAABAuUEBMBybloRIECAAAECBLIVEACzLZ2OEyBAgAABAgTKCQiA5dy0IkCAAAECBAhkKyAAZls6HSdAgAABAgQIlBMQAMu5aUWAAAECBAgQyFZAAMy2dDpOgAABAgQIECgnIACWc9OKAAECBAgQIJCtgACYbel0nAABAgQIECBQTkAALOemFQECBAgQIEAgWwEBMNvS6TgBAgQIECBAoJyAAFjOTSsCBAgQIECAQLYCAmC2pdNxAgQIECBAgEA5AQGwnJtWBAgQIECAAIFsBQTAbEun4wQIECBAgACBcgICYDk3rQgQIECAAAEC2QoIgNmWTscJECBAgAABAuUEBMBybloRIECAAAECBLIVEACzLZ2OEyBAgAABAgTKCQiA5dy0IkCAAAECBAhkKyAAZls6HSdAgAABAgQIlBMQAMu5aUWAAAECBAgQyFZAAMy2dDpOgAABAgQIECgnIACWc9OKAAECBAgQIJCtgACYbel0nAABAgQIECBQTkAALOemFQECBAgQIEAgWwEBMNvS6TgBAgQIECBAoJyAAFjOTSsCBAgQIECAQLYCAmC2pdNxAgQIECBAgEA5AQGwnJtWBAgQIECAAIFsBQTAbEun4wQIECBAgACBcgICYDk3rQgQIECAAAEC2QoIgNmWTscJECBAgAABAuUEBMBybloRIECAAAECBLIVEACzLZ2OEyBAgAABAgTKCQiA5dy0IkCAAAECBAhkKyAAZls6HSdAgAABAgQIlBPo+AD4wgsvxIUXXhjDhg2LI488MmbPnh0vv/xyn3u7a9euuOyyy2LEiBFxxBFHxHnnnRdbt26tt/n+978fM2bMiLFjx8Zhhx0WJ510Utx44431+UYIECBAgAABAv1ZoOMDYAp/jz/+eKxatSruu++++M53vhOXXnppnzX5/Oc/H9/4xjdi+fLl8e1vfzs2b94cn/jEJ+pt1q1bF8ccc0z83d/9XbHuq6++OubPnx8333xzfRkjBAgQIECAAIH+KjCgUh06deee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/fXBqEUKJuns676P9J+uoadAujdi+sUm3Q8wnflI98asDel2RjNnzqy9LJ7/4R/+ofKrv/qrxfLpzNK//uu/9pifbgVzzTXXVKoXghThsvoRfOXpp5/usYwXvQu0sha1/8P2/RlIrxv/X+u9J6a2sha9aQqAvam0f9qA1IXqD4mBAAECBAgQIECgSwRcBdwlhbabBAgQIECAAIGagABYk/BMgAABAgQIEOgSAQGwSwptNwkQIECAAAECNQEBsCbhmQABAgQIECDQJQICYJcU2m4SIECAAAECBGoCAmBNwjMBAgQIECBAoEsEBMAuKbTdJECAAAECBAjUBATAmoRnAgQIECBAgECXCAiAXVJou0mAAAECBAgQqAkIgDUJzwQIECBAgACBLhEQALuk0HaTAAECBAgQIFATEABrEp4JECBAgAABAl0iIAB2SaHtJgECBAgQIECgJiAA1iQ8EyBAgAABAgS6REAA7JJC200CBAgQIECAQE1AAKxJeCZAgAABAgQIdImAANglhbabBAgQIECAAIGagABYk/BMgAABAgQIEOgSAQGwSwptNwkQIECAAAECNQEBsCbhmQABAgQIECDQJQICYJcU2m4SIECAAAECBGoCAmBNwjMBAgQIECBAoEsEBMAuKbTdJECAAAECBAjUBATAmoRnAgQIECBAgECXCAiAXVJou0mAAAECBAgQqAkIgDUJzwQIECBAgACBLhEQALuk0HaTAAECBAgQIFATEABrEp4JECBAgAABAl0iIAB2SaHtJgECBAgQIECgJiAA1iQ8EyBAgAABAgS6REAA7JJC200CBAgQIECAQE3g/wO/ORPDo3iMMQAAAABJRU5ErkJggg=="
            ),
            plotImageUrl: null,
            fileTypeToPlot: null,
            xAxisIndex: 0,
            yAxisIndeces: [],
            flags: "",
            plotParameter: "new",
            plotParameters: [],
            plottedFiles: [],
            plottableFiles: {},
            plottableFilenames: []
        }
    },
    computed: {
        sortedPlotFiles() {
            if (!(this.fileTypeToPlot in this.plottableFiles))
                return [];

            let plotFiles = this.plottableFiles[this.fileTypeToPlot];

            plotFiles.sort((a, b) => a.name.localeCompare(b.name));

            return plotFiles;
        },
        plottablePrefixNames() {
            let plottableRunNames = {};

            if (this.fileTypeToPlot == null)
                return null;

            for (const plottable of this.plottableFiles[this.fileTypeToPlot]) {
                plottableRunNames[plottable.prefix] = plottable.runName;
            }
            
            return Object.entries(plottableRunNames)
        },
        plotObject() {
            return plotObjects[this.fileTypeToPlot]
        },
        file() {
            return this.$store.getters["input/getFile"];
        },
        xAxisOptions() {
            if (this.plotObject != null)
                return this.plotObject.xOptions;
            else 
                return []
        },
        yAxisOptions() {
            if (this.plotObject != null)
                return this.plotObject.yOptions;
            else 
                return []
        },
        canCreate() {
            return this.yAxisIndeces.length > 0 && this.yAxisIndeces.indexOf(this.xAxisIndex) == -1 && this.plottedFiles.length > 0
        }
    },
    methods: {
        async fetchOutputs() {
            let plottableFiles = {};

            for (const run of this.selectedRuns) {
                for (const input of run.inputs) {
                    let properFilename = InputService.getProperFilenameOf(input.name, Object.keys(plotObjects));

                    if (properFilename != null) {
                        plottableFiles[properFilename] = plottableFiles[properFilename] || [];
                        plottableFiles[properFilename].push({
                            runName: run.name,
                            name: `Run${run.id}__${outputFile.name}`, 
                            dataProductURI: outputFile.dataProductURI,
                            prefix: `Run${run.id}`
                        });
                    }
                }

                const outputs = await InputService.fetchOutputs(run.id);

                for (const outputFile of outputs) {
                    let properFilename = InputService.getProperFilenameOf(outputFile.name, Object.keys(plotObjects));

                    if (properFilename != null) {
                        plottableFiles[properFilename] = plottableFiles[properFilename] || [];
                        plottableFiles[properFilename].push({
                            runName: run.name,
                            name: `Run${run.id}__${outputFile.name}`, 
                            dataProductURI: outputFile["data-product-uri"],
                            prefix: `Run${run.id}`
                        });
                    }
                }
            }

            this.plottableFiles = plottableFiles;
            this.plottableFilenames = Object.keys(this.plottableFiles);

            if (!(this.fileTypeToPlot in this.plottableFiles))
                this.fileTypeToPlot = Object.keys(this.plottableFiles)[0];

            await this.fetchParameters();
        },
        async fetchParameters() {
            this.$store.commit("loading/START", { key: "plot", message: "Fetching Plot Parameters" });
            
            try {
                this.plotParameters = (await this.$store.dispatch("plotParameters/fetchPlotParameters"))
                    .map(({xAxis, yAxes, flags}) => {
                        const obj = this.plotObject.xOptions.find(({ value }) => value == xAxis);
                        const xAxisName = (obj != null && "text" in obj) ? obj.text : "--";

                        const yAxisNames = yAxes.split(",").map(yAxis => {
                            const obj = this.plotObject.yOptions.find(({ value }) => value == parseInt(yAxis));

                            return (obj != null && "text" in obj) ? obj.text : "--";
                        });

                        return {
                            text: `x=${xAxisName}, y=[${yAxisNames}]${(flags) ? ", flags=" : ""}${flags}`,
                            value: `x=${xAxis}, y=[${yAxes}]${(flags) ? ", flags=" : ""}${flags}`
                        };
                    });
            } catch (error) {
                eventBus.$emit("error", { name: "Error while trying to load plot parameters", error });
            }

            this.$store.commit("loading/STOP", { key: "plot", message: "Fetching Plot Parameters" });
        },
        async fetchPlot() {
            this.$store.commit("loading/START", { key: "plot", message: "Creating plot" });
            PlotService.plotSelectedRuns({
                plotfiles: this.plottedFiles,
                xAxis: "" + this.xAxisIndex,
                yAxis: this.yAxisIndeces.join(","),
                flags: this.flags
            }).then(({plotImageUrl, output, userGuidance}) => {
                if (plotImageUrl) {
                    this.plotImageUrl = plotImageUrl;
                } else {
                    const error = new Error();
                    error.stack = `[Output]\n${output}`;
                    throw error;
                }

                this.$store.commit("loading/STOP", { key: "plot", message: "Creating plot" });
            }).catch(error => {
                eventBus.$emit("error", { name: "Error while creating the plot", error });
            })
        }
    },
    watch: { 
        plotParameter() {
            if (this.plotParameter != "new") {
                this.xAxisIndex = parseInt(this.plotParameter.match(/x=([0-9]+)/)[1]);
                this.yAxisIndeces = this.plotParameter.match(/y=\[([0-9,]+)\]/)[1].split(",").map(y => parseInt(y));
                this.flags = (this.plotParameter.match(/flags=([^]+)/) || ["", ""])[1];
            }
        },
        selectedRuns() {
            this.fetchOutputs();
        }
    },
    mounted() {
        this.fetchOutputs();
    }
}
</script>
