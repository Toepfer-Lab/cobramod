import {galata, test} from "@jupyterlab/galata";
import {expect} from "@playwright/test";
import * as path from 'path';


test.describe("Visual regression for 'cobramod.escher'.", () => {
  test.beforeEach(async ({ page, tmpPath }) => {
    await page.contents.uploadDirectory(
      path.resolve(__dirname, "./notebooks"),
      tmpPath
    );
    await page.filebrowser.openDirectory(tmpPath);
  });

  test('Using escher integration and compare with reference.', async ({
    page,
    tmpPath,
  }) => {
    const notebook = 'cobramod-escher.ipynb';
    await page.notebook.openByPath(`${tmpPath}/${notebook}`);
    await page.notebook.activate(notebook);

    const captures = new Array<Buffer>();
    let cellCount = 0

    await page.notebook.runCellByCell({
      onAfterCellRun: async (cellIndex: number) => {
        const cell = await page.notebook.getCellOutput(cellIndex);
        if (cell) {
          await new Promise(f => setTimeout(f, 1000));
          captures.push(await cell.screenshot());
          cellCount++
        }
      },
    });

    await page.notebook.save();

    for (let i = 0; i < cellCount; i++) {
      const image = `widgets-cell-${i}.png`;
      expect.soft(captures[i]).toMatchSnapshot(image);
    }
  });
});